#!/usr/bin/env python
from distutils.core import setup, Extension
import numpy

setup( name = "xrr", 
       version = "0.9.07", 
       ext_modules = [Extension("pyxrr.xrr", ["pyxrr/xrr_sc.c"], include_dirs = [numpy.get_include()])],
       packages = ["pyxrr"],
       package_data={'pyxrr': ['materials.sqlite']},
       author = "Carsten Richter", 
       author_email = "carsten.richter@desy.de",
       description = "Contains only one function, which is for calculating X-Ray-Reflectivity (xrr)",
       long_description = "Contains only one function, which is for calculating X-Ray-Reflectivity (xrr).\n\
                           It is supposed to be used by the python xrr wrapper (pyxrr)\n\n\
                           See xrr.reflectivity for more details.",
       )