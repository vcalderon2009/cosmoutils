#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import glob
import os
import sys
from setuptools.command.test import test as TestCommand
import sys

# Enforce Python version check - this is the same check as in __init__.py but
# this one has to happen before importing ah_bootstrap.
if sys.version_info < tuple((int(val) for val in "3.5".split('.'))):
    sys.stderr.write("ERROR: cosmo_utils requires Python {} or later\n".format(3.5))
    sys.exit(1)

from setuptools import setup, find_packages

# Get some values from the setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

conf = ConfigParser()
conf.read(['setup.cfg'])
metadata = dict(conf.items('metadata'))

PACKAGENAME = metadata.get('package_name', 'cosmo_utils')
DESCRIPTION = metadata.get('description', 'Repository with scripts used in my LSS research')
AUTHOR = metadata.get('author', 'Victor Calderon')
AUTHOR_EMAIL = metadata.get('author_email', '')
LICENSE = metadata.get('license', 'unknown')
URL = metadata.get('url', 'https://github.com/vcalderon2009/cosmo_utils')

# order of priority for long_description:
#   (1) set in setup.cfg,
#   (2) load LONG_DESCRIPTION.rst,
#   (3) load README.rst,
#   (4) package docstring
readme_glob = 'README*'
_cfg_long_description = metadata.get('long_description', '')
if _cfg_long_description:
    LONG_DESCRIPTION = _cfg_long_description

elif os.path.exists('LONG_DESCRIPTION.rst'):
    with open('LONG_DESCRIPTION.rst') as f:
        LONG_DESCRIPTION = f.read()

elif len(glob.glob(readme_glob)) > 0:
    with open(glob.glob(readme_glob)[0]) as f:
        LONG_DESCRIPTION = f.read()

else:
    # Get the long description from the package's docstring
    __import__(PACKAGENAME)
    package = sys.modules[PACKAGENAME]
    LONG_DESCRIPTION = package.__doc__

# VERSION should be PEP440 compatible (http://www.python.org/dev/peps/pep-0440)
VERSION = metadata.get('version', '0.0.dev0')

# Indicates if this version is a release version
RELEASE = 'dev' not in VERSION


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

# Compiling Cython
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools.extension import Extension
from find_cython import get_extensions

# extensions = get_extensions()

## Main Setup for the Package

setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      url=URL,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      packages=find_packages(),
      install_requires=[s.strip() for s in metadata.get('install_requires', 'astropy').split(',')],
      setup_requires=[s.strip() for s in metadata.get('setup_requirements', 'astropy').split(',')],
      tests_require=[s.strip() for s in metadata.get('test_requirements', 'astropy').split(',')],
      test_suite='tests',
      cmdclass = {'test': PyTest, 'build_ext': build_ext},
      zip_safe=False,
      # ext_modules = cythonize(extensions)
      )
