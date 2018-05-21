#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['get_extensions']

import os
import glob
import numpy as np
from distutils.extension import Extension

try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

def get_extensions():
    """
    Gets Cython Extensions for given Cython files.

    Returns
    --------
    extension_arr : list
        List of extensions for the Cython scripts.
    """
    # Get some values from the setup.cfg
    conf = ConfigParser()
    conf.read(['setup.cfg'])
    metadata = dict(conf.items('metadata'))
    pkgname  = metadata.get('package_name')
    ## Finding Cython modules
    cython_arr = glob.glob('{0}/**/*.pyx'.format(pkgname),
        recursive=True)
    # Extra args
    include_dirs       = [np.get_include()]
    libraries          = []
    language           = 'c++'
    extra_compile_args = ['-Ofast']
    ##
    ## Looping over Cython files
    extension_arr = []
    for file in cython_arr:
        PATH_TO_PKG   = os.path.relpath(os.path.dirname(file))
        SOURCES       = os.path.basename(file)
        THIS_PKG_NAME = PATH_TO_PKG.replace('/', '.')
        name = PATH_TO_PKG.replace('/','.') + '.' + SOURCES.replace('.pyx', '')
        source = os.path.join(PATH_TO_PKG, SOURCES)
        # Appending to `extension_arr`
        extension_arr.append(Extension(name=name,
            sources=[source],
            include_dirs=include_dirs,
            libraries=libraries,
            language=language,
            extra_compile_args=extra_compile_args))

    return extension_arr
