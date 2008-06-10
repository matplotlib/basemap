To build docs ....
* install sphinx from svn.
* make _static and sphinxext symlinks here that point to the corresponding
  directories in the matplotlib doc tree, e.g.

  ln -s ../../../matplotlib/doc/_static _static
  ln -s ../../../matplotlib/doc/sphinxext sphinxext

* run make.py.  Entry point is build/html/index.html.
