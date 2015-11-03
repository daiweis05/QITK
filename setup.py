
import sys
try:
    import py2exe
except:
    raw_input('Please install py2exe first...')
    sys.exit(-1)
 
from distutils.core import setup
import shutil
 
sys.argv.append('py2exe')

setup( windows=['QITKV0p7.pyw'])