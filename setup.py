#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, Extension, find_packages
from src.lib.aeroframe_2.__version__ import __version__
import os
import platform


def make_path_os_indep(unix_path):
    """
    Return a file path suitable for current OS

    Args:
        :unix_path: Unix-style file path (string)

    Returns:
        :os_path: OS-independent file path (string)
    """

    os_path = os.path.join(*unix_path.split('/'))
    return os_path


# Numpy must be installed
try:
    import numpy as np
except ImportError:
    msg = f"""
    {'='*80}
    Please install 'numpy' first. Try:

        pip install numpy
    {'='*80}
    """
    raise ImportError(msg)

# # Pytornado must be installed
# try:
#     import pytornado as pyt
# except ImportError:
#     msg = f"""
#     {'='*80}
#     Please install 'Pytornado' first. Try:

#         pip install Pytornado
#     {'='*80}
#     """
#     raise ImportError(msg)

# TODO: add test of FramAT, for Sphepy, Fenics, SU2

# See also: https://github.com/kennethreitz/setup.py/blob/master/setup.py
###
NAME = 'aeroframe_2'
VERSION = __version__
AUTHOR = 'Jean-Philippe Kuntzer'
EMAIL = 'jean-philippe.kuntzer19@imperial.ac.uk'
DESCRIPTION = 'Fluid structure interaction module for SU2 and pytornado'
URL = 'https://github.com/jphkun/aeroframe_2'
REQUIRES_PYTHON = '>=3.6.0'
REQUIRED = [
    'airfoils>=0.2.0',
    'aeroframe>=0.0.4',
    'ambiance',
    'numpy',
    'scipy',
    'matplotlib>=3.0.2',
    'commonlibs>=0.5.0',
    'schemadict>=0.0.8',
    # TODO add framat, fenics, SU2, pytornado
]
HERE = os.path.dirname(__file__)
README = os.path.join(HERE, 'README.rst')
PACKAGE_DIR = make_path_os_indep('src/lib/')
LICENSE = 'Apache License 2.0'

# ===== Main executable file =====
SCRIPTS = [make_path_os_indep('src/bin/_aeroframe_2_exe.py')]

# Windows
if platform.system().lower() == 'windows':
    # Use BAT file as wrapper, see file header for reason
    # TODO add specific file
    SCRIPTS.append(make_path_os_indep('src/bin/aeroframe.bat'))
# Linux and MacOs
else:
    SCRIPTS.append(make_path_os_indep('src/bin/aeroframe_2'))

# Note: When downloading https://github.com/airinnova/pytornado/archive/master.zip from Windows,
# the README file is not included for some reason; but we do not need to throw error
if os.path.isfile(README):
    with open(README, "r") as fp:
        long_description = fp.read()
else:
    print("WARNING: Could not find README'")
    long_description = ''

# # Extension modules
# # See
# # * https://docs.python.org/3.7/extending/building.html
# # * https://docs.python.org/3/distutils/setupscript.html
# module1_path_rel = make_path_os_indep('src/lib/pytornado/aero/')
# module1 = Extension(
#         'pytornado.aero.c_vlm',
#         sources=[
#             os.path.join(module1_path_rel, 'c_vlm.cpp'),
#             os.path.join(module1_path_rel, 'c_boundary.cpp'),
#             os.path.join(module1_path_rel, 'c_downwash.cpp'),
#             os.path.join(module1_path_rel, 'c_lattice.cpp'),
#             os.path.join(module1_path_rel, 'c_results.cpp'),
#             ],
#         include_dirs=[
#             module1_path_rel,
#             make_path_os_indep('/usr/include/python3.6'),
#             # Numpy header
#             # '/usr/lib/python3/dist-packages/numpy/core/include/'
#             os.path.join(np.get_include(), 'numpy'),
#             np.get_include(),
#             ]
#         )

setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    long_description=long_description,
    url=URL,
    # ext_modules=[module1],
    include_package_data=True,
    scripts=SCRIPTS,
    package_dir={'': PACKAGE_DIR},
    license=LICENSE,
    # packages=[NAME],
    packages=find_packages(where=PACKAGE_DIR),
    python_requires=REQUIRES_PYTHON,
    install_requires=REQUIRED,
    # See: https://pypi.org/classifiers/
    classifiers=[
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.6',
        "License :: OSI Approved :: Apache Software License",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
