import sys, os
from cx_Freeze import setup, Executable

os.environ['TCL_LIBRARY'] = r'c:\Users\Nick.Forfinski-Sarko\AppData\Local\Continuum\anaconda3\envs\cblue_py36\tcl\tcl8.6'
os.environ['TK_LIBRARY'] = r'c:\Users\Nick.Forfinski-Sarko\AppData\Local\Continuum\anaconda3\envs\cblue_py36\tcl\tk8.6'

__version__ = 'v2.1.0-rc1'

if sys.platform == "win32":
    base = "Win32GUI"

include_files = ['cBLUE_ASCII.txt',
                 'cBLUE_ASCII_finished.txt',
                 'cblue_configuration.json',
                 'cBLUE_icon.ico',
                 'cBLUE_readme.gif',
                 'cBLUE_splash.gif',
                 'lookup_tables']

dll_dir = r'..\DLLs'
for dll in os.listdir(dll_dir):
    include_files.append(os.path.join(dll_dir, dll))

excludes = []
includes = ['numpy.core._methods']
packages = []

setup(
    name = 'cBLUE',
    description='NOAA RSD TPU Tool',
    version=__version__,
    options = {'build_exe': {
    'packages': packages,
    'include_files': include_files,
    'includes': includes,
    'excludes': excludes,
    'include_msvcr': True,
    }},
    executables = [Executable('cBLUEApp.py')]
    )