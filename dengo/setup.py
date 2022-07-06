#!/usr/bin/env python
import os
import sys

import setuptools


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration("dengo", parent_package, top_path)
    config.add_data_dir("solvers")
    config.add_data_dir("templates")
    return config
