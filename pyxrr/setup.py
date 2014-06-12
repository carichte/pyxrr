#!/usr/bin/env python
from distutils.core import setup, Extension
import numpy
import sys

if len(sys.argv)<2:
    print("see install.txt for installation instructions.")

if "--single-core" in sys.argv[1:]:
    ext_modules = [Extension("pyxrr.xrr", ["pyxrr/xrr_sc.c"], include_dirs = [numpy.get_include()])]
    sys.argv.pop(sys.argv.index("--single-core"))
else:
    ext_modules = [Extension("pyxrr.xrr",
                             ["pyxrr/xrr.c"], 
                             include_dirs = [numpy.get_include()], 
                             extra_compile_args = ["-fopenmp"], 
                             extra_link_args = ["-fopenmp"])]

setup( name = "xrr", 
       version = "0.9.08",
       ext_modules = ext_modules,
       packages = ["pyxrr"],
       package_data={'pyxrr': ['materials.sqlite',
                               'locale/en/LC_MESSAGES/*',
                               'locale/de/LC_MESSAGES/*'
                               ]},
       author = "Carsten Richter", 
       author_email = "carsten.richter@desy.de",
       description = "Contains only one function, which is for calculating X-Ray-Reflectivity (xrr)",
       long_description = "Contains only one function, which is for calculating X-Ray-Reflectivity (xrr).\n\
                           It is supposed to be used by the python xrr wrapper (pyxrr)\n\n\
                           See xrr.reflectivity for more details.",
     )

