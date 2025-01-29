#!/usr/bin/env python

import os
from setuptools import setup

setup(
    name = 'spp-dcj',
    version = '2.0',
    url = 'https://github.com/marschall-lab/spp_dcj_v2',
    author = 'Leonard Bohnenkämper, Daniel Dörr',
    author_email = 'lbohnenk@CeBiTec.Uni-Bielefeld.DE, daniel.doerr@HHU.de',
    description = 'A tool for solving the small parsimony problem for natural genomes',
    long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
    license = 'MIT',
    keywords = 'double-cut-and-join genome rearrangement bioinformatics',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=['spp_dcj'],
    scripts=['scripts/spp-dcj'],
)
