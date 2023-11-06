from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from ampli_clean import __version__, _program

setup(
    name='ampli_clean',
    version='0.0.1',    
    description='ampli_clean',
    url='https://github.com/Desperate-Dan/ampli_clean',
    author='Daniel Maloney',
    author_email='dmaloney@ed.ac.uk',
    packages=find_packages(),
    scripts=['ampli_clean/ampli_clean.py'],
    install_requires=['biopython',
                        "pysam",
                        "numpy",
                        "mappy"
                      ],
    entry_points="""
    [console_scripts]
    {program} = ampli_clean:main
    """.format(program = _program)
)
