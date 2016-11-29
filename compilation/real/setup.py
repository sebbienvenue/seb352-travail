#!/usr/bin/env python

"""
setup.py file for SWIG sobol
"""

from distutils.core import setup, Extension


sobol_module = Extension('_sobol',
                           sources=['sobol_wrap.cxx', 'sobol.cpp'],
                           )

setup (name = 'sobol',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [sobol_module],
       py_modules = ["sobol"],
       )
