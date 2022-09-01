import os
from setuptools import setup, find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
# def read(fname):
#     return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'pyaoa',
    version = '0.1',  # Ideally should be same as your GitHub release tag varsion
    description = 'Python package to automatically analyze objects in external aerodynamic flow for a range of angle-of-attack',
    author = "inflowencer",
    author_email= "jason7bernhard@gmail.com",
    keywords = "cfd pyaoa aerodynamics",
    url = "http://packages.python.org/pyaoa",
    # download_url = 'https://github.com/inflowencer/pyaoa/archive/refs/tags/v0.1.tar.gz',
    install_requires=[
        "rich",
        "numpy",
        "pandas",
        "matplotlib",
        "pyyaml"
    ],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    # long_description=read('README.md'),
    # classifiers=[
    #     "Development Status :: 3 - Alpha",
    #     "Topic :: Utilities",
    #     # "License :: OSI Approved :: BSD License",
    # ],
)