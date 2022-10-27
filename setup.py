# -*- coding: utf-8 -*-
"""Setup file for qihmc
"""
__author__ = "None"

from os import path
from distutils.extension import Extension

from setuptools import setup, find_packages

from Cython.Build import cythonize
from numpy import get_include

ROOT = path.abspath(path.dirname(__file__))
SRC = path.join(ROOT, "src")

SOURCES = [
    path.join(SRC, "quboIsingHMC.cpp"),
    path.join(SRC, "_qubo_ising_hmc.pyx"),
]


EXTENSIONS = Extension(
    name="qihmc.qihmc_cc",
    sources=SOURCES,
    include_dirs=[SRC, get_include()],
    language="c++14",
    extra_compile_args=["-std=c++11", "-ffast-math", "-O2"],  # "-fopenmp"],
    # extra_link_args=["-fopenmp"],
)

with open(path.join(ROOT, "README.md"), encoding="utf-8") as inp:
    LONG_DESCRIPTION = inp.read()

setup(
    description=None,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url=None,
    author=__author__,
    author_email=None,
    keywords=[],
    packages=find_packages(exclude=["docs", "tests"]),
    package_dir={"": "src"},
    ext_modules=cythonize([EXTENSIONS]),
)
