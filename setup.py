import os, os.path
import sys
import time
import subprocess

import setuptools

VERSION="0.1"

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    config.add_subpackage('dengo','dengo')

    return config

def setup_package():

    from numpy.distutils.core import setup

    setup(
        name = "dengo",
        version = VERSION,
        description = "A framework for creating chemical and radiative cooling"
                    + "reaction networks",
        classifiers = [ "Development Status :: 2 - Pre-Alpha",
                        "Environment :: Console",
                        "Intended Audience :: Science/Research",
                        "License :: OSI Approved :: GNU General Public License (GPL)",
                        "Operating System :: MacOS :: MacOS X",
                        "Operating System :: POSIX :: Linux",
                        "Programming Language :: C",
                        "Programming Language :: Python",
                        "Topic :: Scientific/Engineering :: Astronomy",
                        "Topic :: Scientific/Engineering :: Physics", ],
        keywords='astronomy astrophysics chemistry',
        author="Matthew J. Turk",
        author_email="matthewturk@gmail.com",
        url = "http://bitbucket.org/MatthewTurk/dengo",
        license="GPL-3",
        configuration=configuration,
        zip_safe=False,
    )
    return

if __name__ == '__main__':
    setup_package()
