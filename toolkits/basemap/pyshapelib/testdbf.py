# Copyright (C) 2003 by Intevation GmbH
# Authors:
# Bernhard Herzog <bh@intevation.de>
#
# This program is free software under the LGPL (>=v2)
# Read the file COPYING coming with the software for details.

"""Test cases for the dbflib python bindings"""

__version__ = "$Revision$"
# $Source$
# $Id$

import unittest
import dbflib

class TestDBF(unittest.TestCase):

    def test_add_field(self):
        """Test whethe add_field reports exceptions"""
        dbf = dbflib.create("test.dbf")
        # For strings the precision parameter must be 0
        self.assertRaises(RuntimeError,
                          dbf.add_field, "str", dbflib.FTString, 10, 5)


if __name__ == "__main__":
    unittest.main()
