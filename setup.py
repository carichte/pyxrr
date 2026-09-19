#!/usr/bin/env python
import os, sys
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext



if len(sys.argv)<2:
    print("see install.txt for installation instructions.")


class build_shared_lib(build_ext):
    """Compile libxrr.c as a plain shared library (loaded via ctypes)."""
    def build_extension(self, ext):
        target = os.path.join(self.build_lib, 'pyxrr', 'libxrr')
        os.makedirs(os.path.dirname(target), exist_ok=True)

        if self.compiler.compiler_type == 'msvc':
            compile_args = ['/openmp']
            link_args = []
        else:
            compile_args = ext.extra_compile_args
            link_args = list(ext.extra_link_args)
            if sys.platform == 'win32':
                link_args += ['-static-libgcc',
                              '-Wl,-Bstatic',
                              '-lgomp', '-lwinpthread', '-ldl',
                              '-Wl,-Bdynamic']
        
        objects = self.compiler.compile(
            ext.sources,
            output_dir=self.build_temp,
            extra_postargs=compile_args,
            include_dirs=ext.include_dirs)

        if sys.platform == 'win32':
            target += '.dll'
        else:
            target += '.so'

        self.compiler.link_shared_object(
            objects, target,
            extra_postargs=link_args)


ext_modules = [Extension("pyxrr.libxrr",
               ["pyxrr/libxrr.c"],
               include_dirs=[],
               extra_compile_args=["-fopenmp"],
               extra_link_args=["-fopenmp"])]


setup(
    name = "pyxrr", 
    version = "1.0.0",
    ext_modules = ext_modules,
    cmdclass={'build_ext': build_shared_lib},
    packages = ["pyxrr"],
    package_data={'pyxrr': [
        'f1f2.h5',
        'locale/en/LC_MESSAGES/*',
        'locale/de/LC_MESSAGES/*'
    ]},
    author = "Carsten Richter", 
    author_email = "carsten.richter@ikz-berlin.de",
    description = "Contains a function for calculating X-Ray-Reflectivity (libxrr)",
    long_description = "Contains a function for calculating X-Ray-Reflectivity (libxrr).\n\
    It is supposed to be used by the python xrr wrapper (pyxrr).",
    install_requires=[
        'numpy',
        'lmfit',
        'appdirs',
        'matplotlib',
        'scipy',
        'h5py',
        'pandas',
    ],
)

