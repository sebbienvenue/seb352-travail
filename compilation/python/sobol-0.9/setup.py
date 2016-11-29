#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name='sobol',
      version='0.9',
      author='Victor Ng',
      author_email='victor@crankycoder.com',
      description='SOBOL quasi random number sequence generator',
      url="http://github.com/crankycoder/sobol_seq  ",
      license='MIT License',
      packages=find_packages(exclude=['tests']))
