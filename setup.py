# -*- coding: utf-8 -*-
"""Setup file for qubo_ising_hmc
"""
__author__ = "None"
__version__ = "0.1.0"

from os import path
from distutils.extension import Extension

from setuptools import setup, find_packages

from Cython.Build import cythonize
from numpy import get_include

ROOT = path.abspath(path.dirname(__file__))
SRC = path.join(ROOT, "src")
PY_SRC = path.join(ROOT, "qubo_ising_hmc")

SOURCES = [
    path.join(SRC, "quboIsingHMC.cpp"),
    path.join(PY_SRC, "_qubo_ising_hmc.pyx"),
]

EXTENSIONS = Extension(
    name="qubo_ising_hmc._qubo_ising_hmc",
    sources=SOURCES,
    include_dirs=[SRC, get_include()],
    language="c++",
)

CWD = path.abspath(path.dirname(__file__))

with open(path.join(CWD, "README.md"), encoding="utf-8") as inp:
    LONG_DESCRIPTION = inp.read()

with open(path.join(CWD, "requirements.txt"), encoding="utf-8") as inp:
    REQUIREMENTS = [el.strip() for el in inp.read().split("\n")]

setup(
    name="qubo_ising_hmc",
    version=__version__,
    description=None,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url=None,
    author=__author__,
    author_email=None,
    keywords=[],
    packages=find_packages(exclude=["docs", "tests"]),
    install_requires=REQUIREMENTS,
    ext_modules=cythonize([EXTENSIONS]),
)
