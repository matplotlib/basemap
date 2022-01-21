To build docs ....

* install sphinx.

* run 'python make.py html'.  Entry point is build/html/index.html.

* to update on github site:
  checkout gh-pages basemap branch in ../../basemap.git.web
  (cd ../..; git clone git@github.com:matplotlib/basemap.git -b gh-pages basemap.git.web)
  rsync -vaz build/html/ ../../basemap.git.web
  go to basemap.git.web, commit and push the updates.
