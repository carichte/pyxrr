#!/usr/bin/env python
from distutils.core import setup, Extension
import numpy
import sys



#raise Exception("You downloaded the (not working/under construction) developement version. "
#                "Please refer to http://carichte.github.io/pyxrr/ to download "
#                "the stable Version.")


if len(sys.argv)<2:
    print("see install.txt for installation instructions.")

ext_modules = [Extension("pyxrr.libxrr",
               ["pyxrr/libxrr.c"],
               include_dirs = [numpy.get_include()], 
               extra_compile_args = ["-fopenmp"],
               extra_link_args = ["-fopenmp"])]


setup( name = "pyxrr", 
       version = "1.0.0",
       ext_modules = ext_modules,
       packages = ["pyxrr"],
       package_data={'pyxrr': ['f1f2.h5',
                               'locale/en/LC_MESSAGES/*',
                               'locale/de/LC_MESSAGES/*'
                               ]},
       author = "Carsten Richter", 
       author_email = "carsten.richter@physik.tu-freiberg.de",
       description = "Contains a function for calculating X-Ray-Reflectivity (libxrr)",
       long_description = "Contains a function for calculating X-Ray-Reflectivity (libxrr).\n\
                           It is supposed to be used by the python xrr wrapper (pyxrr).",
     )

