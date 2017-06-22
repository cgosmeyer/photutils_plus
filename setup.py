#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

setup(name = 'photutils_plus',
      description = 'Photutils extension package',
      author = 'C.M. Gosmeyer',
      url = 'https://github.com/cgosmeyer/photutils_plus',
      packages = find_packages(),
      install_requires = ['astropy', 'numpy', 'matplotlib', 'photutils']
     )