import sys
import inspect
import compiler
from urllib import quote

from dap.lib import encode_atom


class ASTVisitor(compiler.visitor.ASTVisitor):
    def __init__(self, scope, seq):
        compiler.visitor.ASTVisitor.__init__(self)
        self.scope = scope
        self.seq = seq

        self.filters = []

    def visitGenExprFor(self, node):
        # Check if we're iterating over the sequence.
        if self.eval(node.iter) == self.seq:
            # Add sequence with assigned name to the scope.
            self.scope[node.assign.name] = self.seq

        # Carry on.
        self.visit(node.assign)
        self.visit(node.iter)
        for if_ in node.ifs: 
            self.visit(if_)

    def visitListCompFor(self, node):
        # Check if we're iterating over the sequence.
        if self.eval(node.list) == self.seq:
            # Add sequence with assigned name to the scope.
            self.scope[node.assign.name] = self.seq

        # Carry on.
        self.visit(node.assign)
        self.visit(node.list)
        for if_ in node.ifs:
            self.visit(if_)

    def visitCompare(self, node):
        # Convert multiple comparisons like 1 < a < 2 to
        # (1 < a) and (a < 2).
        if len(node.ops) > 1:
            left = node.expr
            out = []
            for op in node.ops:
                out.append(compiler.ast.Compare(left, [op]))
                left = op[1]
            new_node = compiler.ast.And(out)
            self.visit(new_node)

        else:
            # Simple comparison.
            a, op, b = node.getChildren()
            ops = ['<', '>', '==', '!=', '<=', '>=']
            if op in ops:
                a = self.eval(a)
                b = self.eval(b)
            
                # If the objects have an 'id' attribute we use it
                # instead. 
                if hasattr(a, 'id'): a = a.id
                else: a = quote(encode_atom(a))

                if hasattr(b, 'id'): b = b.id
                else: b = quote(encode_atom(b))

                # Build filter.
                if op == '==': op = '='
                filter_ = '%s%s%s' % (a, op, b)
                self.filters.append(filter_)

    def visitOr(self, node):
        raise Exception('OR not supported by the DAP spec!')

    def eval(self, node):
        """
        Eval node.

        This is done by converting the node to bytecode and
        eval()ing the bytecode in the instance scope.
        """
        ast = compiler.ast.Expression(node)
        ast.filename = 'dummy'
        c = compiler.pycodegen.ExpressionCodeGenerator(ast)
        obj = eval(c.getCode(), self.scope)

        return obj


def get_filters(seq):
    # We need to get the stack where the sequence is being iterated.
    # This is called from genexp --> __iter__ --> get_filters, so
    # we go up 2 frames in the stack.
    frame = sys._getframe(2)

    # Inspect frame.
    fname, lineno, func, src, index = inspect.getframeinfo(frame)
    scope = frame.f_globals
    
    # Fix source.
    if src:
        src = src[0].strip()
        if src.endswith(':'): src = '%s pass' % src
    
    # Build and walk AST to parse filters.
    visitor = ASTVisitor(scope, seq)
    try:
        ast = compiler.parse(src)
        compiler.walk(ast, visitor)
    except:
        pass

    return visitor.filters
