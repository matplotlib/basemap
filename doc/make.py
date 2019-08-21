#!/usr/bin/env python
from __future__ import print_function
import fileinput
import glob
import os
import shutil
import sys

def html():
    os.system('sphinx-build -b html -d build/doctrees . build/html')

def latex():
    if sys.platform != 'win32':
        # LaTeX format.
        os.system('sphinx-build -b latex -d build/doctrees . build/latex')

        # Produce pdf.
        os.chdir('build/latex')

        # Copying the makefile produced by sphinx...
        os.system('pdflatex Basemap.tex')
        os.system('pdflatex Basemap.tex')
        os.system('makeindex -s python.ist Basemap.idx')
        os.system('makeindex -s python.ist modBasemap.idx')
        os.system('pdflatex Basemap.tex')

        os.chdir('../..')
    else:
        print('latex build has not been tested on windows')

def clean():
    shutil.rmtree('build')

def all():
    html()
    latex()


funcd = {
         'html':html,
         'latex':latex,
         'clean':clean,
         'all':all,
         }


if len(sys.argv)>1:
    for arg in sys.argv[1:]:
        func = funcd.get(arg)
        if func is None:
            raise SystemExit('Do not know how to handle %s; valid args are'%(
                    arg, list(funcd.keys())))
        func()
else:
    all()
