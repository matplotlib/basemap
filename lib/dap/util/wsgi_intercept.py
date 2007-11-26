"""
wsgi_intercept.WSGI_HTTPConnection is a replacement for
httplib.HTTPConnection that intercepts certain HTTP connections into a
WSGI application.

Use 'add_wsgi_intercept' and 'remove_wsgi_intercept' to control this behavior.
"""

import sys
from httplib import HTTPConnection
from cStringIO import StringIO
import traceback

debuglevel = 0
# 1 basic
# 2 verbose

####

#
# Specify which hosts/ports to target for interception to a given WSGI app.
#
# For simplicity's sake, intercept ENTIRE host/port combinations;
# intercepting only specific URL subtrees gets complicated, because we don't
# have that information in the HTTPConnection.connect() function that does the
# redirection.
#
# format: key=(host, port), value=(create_app, top_url)
#
# (top_url becomes the SCRIPT_NAME)

_wsgi_intercept = {}

def add_wsgi_intercept(host, port, app_create_fn, script_name=''):
    """
    Add a WSGI intercept call for host:port, using the app returned
    by app_create_fn with a SCRIPT_NAME of 'script_name' (default '').
    """
    _wsgi_intercept[(host, port)] = (app_create_fn, script_name)

def remove_wsgi_intercept(host, port):
    """
    Remove the WSGI intercept call for (host, port).
    """
    key = (host, port)
    if _wsgi_intercept.has_key(key):
        del _wsgi_intercept[key]

#
# make_environ: behave like a Web server.  Take in 'input', and behave
# as if you're bound to 'host' and 'port'; build an environment dict
# for the WSGI app.
#
# This is where the magic happens, folks.
#

def make_environ(inp, host, port, script_name):
    """
    Take 'inp' as if it were HTTP-speak being received on host:port,
    and parse it into a WSGI-ok environment dictionary.  Return the
    dictionary.

    Set 'SCRIPT_NAME' from the 'script_name' input, and, if present,
    remove it from the beginning of the PATH_INFO variable.
    """
    #
    # parse the input up to the first blank line (or its end).
    #

    environ = {}
    
    method_line = inp.readline()
    
    content_type = None
    content_length = None
    cookies = []
    
    for line in inp:
        if not line.strip():
            break

        k, v = line.strip().split(':', 1)
        v = v.lstrip()

        #
        # take care of special headers, and for the rest, put them
        # into the environ with HTTP_ in front.
        #

        if k.lower() == 'content-type':
            content_type = v
        elif k.lower() == 'content-length':
            content_length = v
        elif k.lower() == 'cookie' or k.lower() == 'cookie2':
            cookies.append(v)
        else:
            h = k.upper()
            h = h.replace('-', '_')
            environ['HTTP_' + h] = v
            
        if debuglevel >= 2:
            print 'HEADER:', k, v

    #
    # decode the method line
    #

    if debuglevel >= 2:
        print 'METHOD LINE:', method_line
        
    method, url, protocol = method_line.split(' ')

    # clean the script_name off of the url, if it's there.
    if not url.startswith(script_name):
        script_name = ''                # @CTB what to do -- bad URL.  scrap?
    else:
        url = url[len(script_name):]

    url = url.split('?', 1)
    path_info = url[0]
    query_string = ""
    if len(url) == 2:
        query_string = url[1]

    if debuglevel:
        print "method: %s; script_name: %s; path_info: %s; query_string: %s" % (method, script_name, path_info, query_string)

    r = inp.read()
    inp = StringIO(r)

    #
    # fill out our dictionary.
    #
    
    environ.update({ "wsgi.version" : (1,0),
                     "wsgi.url_scheme": "http",
                     "wsgi.input" : inp,           # to read for POSTs
                     "wsgi.errors" : StringIO(),
                     "wsgi.multithread" : 0,
                     "wsgi.multiprocess" : 0,
                     "wsgi.run_once" : 0,
    
                     "REQUEST_METHOD" : method,
                     "SCRIPT_NAME" : script_name,
                     "PATH_INFO" : path_info,
          
                     "SERVER_NAME" : host,
                     "SERVER_PORT" : str(port),
                     "SERVER_PROTOCOL" : protocol,

                     "REMOTE_ADDRESS" : '127.0.0.1',
                     })

    #
    # query_string, content_type & length are optional.
    #

    if query_string:
        environ['QUERY_STRING'] = query_string
        
    if content_type:
        environ['CONTENT_TYPE'] = content_type
        if debuglevel >= 2:
            print 'CONTENT-TYPE:', content_type
    if content_length:
        environ['CONTENT_LENGTH'] = content_length
        if debuglevel >= 2:
            print 'CONTENT-LENGTH:', content_length

    #
    # handle cookies.
    #
    if cookies:
        environ['HTTP_COOKIE'] = "; ".join(cookies)

    if debuglevel:
        print 'WSGI environ dictionary:', environ

    return environ

#
# fake socket for WSGI intercept stuff.
#

class wsgi_fake_socket:
    """
    Handle HTTP traffic and stuff into a WSGI application object instead.

    Note that this class assumes:
    
     1. 'makefile' is called (by the response class) only after all of the
        data has been sent to the socket by the request class;
     2. non-persistent (i.e. non-HTTP/1.1) connections.
    """
    def __init__(self, app, host, port, script_name):
        self.app = app                  # WSGI app object
        self.host = host
        self.port = port
        self.script_name = script_name  # SCRIPT_NAME (app mount point)

        self.inp = StringIO()           # stuff written into this "socket"
        self.results = None             # results from running the app
        self.output = StringIO()        # all output from the app, incl headers

    def makefile(self, *args, **kwargs):
        """
        'makefile' is called by the HTTPResponse class once all of the
        data has been written.  So, in this interceptor class, we need to:
        
          1. build a start_response function that grabs all the headers
             returned by the WSGI app;
          2. create a wsgi.input file object 'inp', containing all of the
             traffic;
          3. build an environment dict out of the traffic in inp;
          4. run the WSGI app & grab the result object;
          5. concatenate & return the result(s) read from the result object.

        @CTB: 'start_response' should return a function that writes
        directly to self.result, too.
        """

        # dynamically construct the start_response function for no good reason.
        def start_response(status, headers, exc_info=None):
            
            # construct the HTTP request.
            self.output.write("HTTP/1.0 " + status + "\n")
            
            for k, v in headers:
                self.output.write('%s: %s\n' % (k, v,))
            self.output.write('\n')

        # construct the wsgi.input file from everything that's been
        # written to this "socket".
        inp = StringIO(self.inp.getvalue())

        # build the environ dictionary.
        environ = make_environ(inp, self.host, self.port, self.script_name)

        # run the application.
        self.result = self.app(environ, start_response)

        ###

        # read all of the results.
        for data in self.result:
            self.output.write(data)

        if debuglevel >= 2:
            print "***", self.output.getvalue(), "***"

        # return the concatenated results.
        return StringIO(self.output.getvalue())

    def sendall(self, str):
        """
        Save all the traffic to self.inp.
        """
        if debuglevel >= 2:
            print ">>>", str, ">>>"

        self.inp.write(str)

    def close(self):
        "Do nothing, for now."
        pass

#
# WSGI_HTTPConnection
#

class WSGI_HTTPConnection(HTTPConnection):
    """
    Intercept all traffic to certain hosts & redirect into a WSGI
    application object.
    """
    def get_app(self, host, port):
        """
        Return the app object for the given (host, port).
        """
        key = (host, int(port))

        app, script_name = None, None
        
        if _wsgi_intercept.has_key(key):
            (app_fn, script_name) = _wsgi_intercept[key]
            app = app_fn()

        return app, script_name        
    
    def connect(self):
        """
        Override the connect() function to intercept calls to certain
        host/ports.
        """
        if debuglevel:
            sys.stderr.write('connect: %s, %s\n' % (self.host, self.port,))
                             
        try:
            (app, script_name) = self.get_app(self.host, self.port)
            if app:
                if debuglevel:
                    sys.stderr.write('INTERCEPTING call to %s:%s\n' % \
                                     (self.host, self.port,))
                self.sock = wsgi_fake_socket(app, self.host, self.port,
                                             script_name)
            else:
                HTTPConnection.connect(self)
                
        except Exception, e:
            if debuglevel:              # intercept & print out tracebacks
                traceback.print_exc()
            raise

### DEBUGGING CODE -- to help me figure out communications stuff. ###

# (ignore me, please)

'''
import socket
    
class file_wrapper:
    def __init__(self, fp):
        self.fp = fp

    def readline(self):
        d = self.fp.readline()
        if debuglevel:
            print 'file_wrapper readline:', d
        return d

    def __iter__(self):
        return self

    def next(self):
        d = self.fp.next()
        if debuglevel:
            print 'file_wrapper next:', d
        return d

    def read(self, *args):
        d = self.fp.read(*args)
        if debuglevel:
            print 'file_wrapper read:', d
        return d

    def close(self):
        if debuglevel:
            print 'file_wrapper close'
        self.fp.close()

class intercept_socket:
    """
    A socket that intercepts everything written to it & read from it.
    """

    def __init__(self):
        for res in socket.getaddrinfo("floating.caltech.edu", 80, 0,
                                      socket.SOCK_STREAM):
            af, socktype, proto, canonname, sa = res
            self.sock = socket.socket(af, socktype, proto)
            self._open = True
            self.sock.connect(sa)
            break

    def makefile(self, *args, **kwargs):
        fp = self.sock.makefile('rb', 0)
        return file_wrapper(fp)

    def sendall(self, str):
        if not self._open:
            raise Exception

        return self.sock.sendall(str)
    
    def close(self):
        self._open = False
'''
