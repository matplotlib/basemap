"""Parsing routines.

This module implements parsers for the dataset descriptors: DAS/DDS, DDX, JSON, etc.
"""

__author__ = "Roberto De Almeida <rob@pydap.org>"

import shlex


class BaseParser(shlex.shlex):
    """A base parser with support functions.
    
    This module is a simple parser based on ``shlex.shlex``. It has 
    support function for peeking and consuming tokens.
    """
    def _consume(self, token):
        """Consume the specified token or raise exception."""
        _token = self.get_token()
        if token.lower() == _token.lower(): return _token

        e = "%s Found '%s' (expected '%s')" % (self.error_leader(self.url), _token, token)
        raise Exception(e)

    def _peek(self):
        """Inspect the next token."""
        _token = self.get_token()
        self.push_token(_token)
        return _token
        
    def _check(self, *tokens):
        """Check for token(s)."""
        _token = self._peek()
        if _token.lower() in [t.lower() for t in tokens]: return _token
        return None
