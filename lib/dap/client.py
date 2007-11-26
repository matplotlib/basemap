__author__ = "Roberto De Almeida <rob@pydap.org>"

import dap.lib
from dap.util.http import openurl
from dap.exceptions import ClientError


def open(url, cache=None, username=None, password=None, verbose=False):
    """Connect to a remote dataset.

    This function opens a dataset stored in a DAP server:

        >>> dataset = open(url, cache=None, username=None, password=None, verbose=False):

    You can specify a cache location (a directory), so that repeated
    accesses to the same URL avoid the network.
    
    The username and password may be necessary if the DAP server requires
    authentication. The 'verbose' option will make pydap print all the 
    URLs that are acessed.
    """
    # Set variables on module namespace.
    dap.lib.VERBOSE = verbose

    if url.startswith('http'):
        for response in [_ddx, _ddsdas]:
            dataset = response(url, cache, username, password)
            if dataset: return dataset
        else:
            raise ClientError("Unable to open dataset.")
    else:
        from dap.plugins.lib import loadhandler
        from dap.helper import walk

        # Open a local file. This is a clever hack. :)
        handler = loadhandler(url)
        dataset = handler._parseconstraints()

        # Unwrap any arrayterators in the dataset.
        for var in walk(dataset):
            try: var.data = var.data._var
            except: pass

    return dataset


def _ddsdas(baseurl, cache, username, password):
    ddsurl, dasurl = '%s.dds' % baseurl, '%s.das' % baseurl

    # Get metadata.
    respdds, dds = openurl(ddsurl, cache, username, password)
    respdas, das = openurl(dasurl, cache, username, password)

    if respdds['status'] == '200' and respdas['status'] == '200':
        from dap.parsers.dds import DDSParser
        from dap.parsers.das import DASParser

        # Build dataset.
        dataset = DDSParser(dds, ddsurl, cache, username, password).parse()

        # Add attributes.
        dataset = DASParser(das, dasurl, dataset).parse()
        return dataset


def _ddx(baseurl, cache, username, password):
    pass
