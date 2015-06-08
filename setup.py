#!/usr/bin/env python

from distutils.core import setup

# Get version number from package
exec(open('illqc/version.py').read())

setup(
    name='illQC',
    version=__version__,
    description='Process shotgun metagenomic DNA data for quality',
    author='Kyle Bittinger',
    author_email='kylebittinger@gmail.com',
    url='https://github.com/PennChopMicrobiomeProgram',
    packages=['illqc'],
    scripts=['scripts/illqc.py']
    )
