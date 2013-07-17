import os, sys
sys.path.append('../snappy/')
import version

spkg_name = 'snappy-' + version.version
os.system('rm -R /tmp/' + spkg_name) 
os.system('hg convert --filemap snappy-sage-filemap ../ /tmp/' + spkg_name)
os.chdir('/tmp/' + spkg_name)
os.system('hg up')
os.system('./fetch_source')
os.chdir('..')
os.system('sage --pkg ' + spkg_name)
raw_input('Hit any key when ready to begin copying to t3m:')
os.system('scp -p ' + spkg_name + '.spkg nmd@shell.math.uic.edu:t3m_web/SnapPy-nest')
