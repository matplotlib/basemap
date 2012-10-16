import glob, os, sys
test_files = glob.glob('*.py')
test_files.remove('run_all.py')
test_files.remove('allskymap.py')
test_files.remove('fcstmaps.py')
test_files.remove('fcstmaps_axesgrid.py')
test_files.remove('testgdal.py')
test_files.remove('animate.py')
test_files.remove('geos_demo_2.py')
test_files.remove('plotsst.py')
test_files.remove('embedding_map_in_wx.py') # requires wx
test_files.remove('plothighsandlows.py') # requires scipy
test_files.remove('lic_demo.py')
test_files.remove('testwmsimage.py')
py_path = os.environ.get('PYTHONPATH')
if py_path is None:
    py_path = '.'
else:
    py_path = os.pathsep.join(['.',py_path])
os.environ['PYTHONPATH'] = py_path

for f in test_files:
    sys.stdout.write( "**********************************************\n")
    ff = os.path.join(sys.path[0],f)
    args = [sys.executable,ff]
    sys.stdout.write("Running %s\n" % f)
    status = os.spawnve(os.P_WAIT,sys.executable,args,os.environ)
    if status:
        sys.stdout.write('TEST FAILURE (status=%s)\n' % (status))
