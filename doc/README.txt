To build docs ....

* install sphinx from svn.

* make _static and sphinxext symlinks here that point to the corresponding
  directories in the matplotlib doc tree, e.g.

  ln -s ../../../matplotlib/doc/_static _static
  ln -s ../../../matplotlib/doc/sphinxext sphinxext

* copy matplotlibrc to users/figures, edit to make default figure size 6,6 
  instead of 8,6.

* run 'python make.py html'.  Entry point is build/html/index.html.

* to update on sf site:
  cd build
  rsync -avz html jswhit,matplotlib@web.sf.net:/home/groups/m/ma/matplotlib/htdocs/basemap/doc/ -essh
