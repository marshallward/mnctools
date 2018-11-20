"""Installation script for mnctools
"""
import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

project_name = 'mnctools'
project_version = __import__(project_name).__version__
project_pkgs = [path for (path, dirs, files) in os.walk(project_name)
                if '__init__.py' in files]
project_readme_fname = 'README.rst'

with open(project_readme_fname) as f:
    project_readme = f.read()

setup(
    name = project_name,
    version = project_version,
    description = 'MNC support module for MITgcm',
    long_description = project_readme,
    author = 'Marshall Ward',
    author_email = 'mnctools@marshallward.org',
    url = 'http://github.com/marshallward/mnctools',

    packages=project_pkgs,
    install_requires=[
        'netCDF4',
        'numpy',
    ],
    classifiers = [
        'Development Status :: 1 - Planning',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Topic :: Utilities',
    ],

    keywords = ['MITgcm'],
)
