#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Module docstring.

This serves as a long usage message.
"""
import sys
import getopt
import os

#http://code.activestate.com/recipes/474083-get-the-path-of-the-currently-executing-python-scr/
# BASE_PATH is the absolute path of ../.. relative to this script location

# add ../../scripts (relative to the file (!) and not to the CWD)
sys.path.append( os.path.join( reduce (lambda l,r: l + os.path.sep + r, os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep )[:-2] ), "share/scripts" ) )

import common

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    # etc., replacing sys.argv with argv in the getopt() call.

    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error, msg:
             raise Usage(msg)
        # more code, unchanged
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

    # print __doc__

    common.printhello()
    common.do()
    cfs = common.CFS()
    print cfs.script_bin_dir
    print cfs.distro_out
#    print cfs.rev
    print cfs.getArchStr()
    print cfs.lib_suffix
    print cfs.exe_bin_dir



if __name__ == "__main__":
    sys.exit(main())
