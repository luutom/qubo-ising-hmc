# -*- coding: utf-8 -*-
"""Setup file for qihmc
"""
__author__ = "None"
__version__ = "0.1.0"

from os import path, environ
import shutil
import subprocess
import pathlib
import platform
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

def findgccbin():
    gcc = environ.get("CC")
    if gcc is None:
        gcc = shutil.which("gcc")
        if gcc is None:
            raise SystemError("Missing gcc, add to path or set env CC")
    print("Found gcc at ", gcc)
    # do test run of gcc
    cmd = gcc + " -v"
    ps = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0].decode("utf-8").split("\n")
    if ps[0] != "Using built-in specs.":
        raise SystemError("Requires gcc. Exiting install")
    p = pathlib.Path(gcc)
    # set some associated env vars for make
    environ["CC"] = p.as_posix()
    environ["CXX"] = p.parents[0].joinpath('g++').as_posix()
    # Fix failure on unbuntu
    # https://stackoverflow.com/questions/45308426/distutils-partially-ignoring-cc-environment-variable-when-building-python-extens
    if platform.system() != 'Darwin':
        environ["LDSHARED"] = environ["CXX"] + " -std=c++14 -pthread -shared"
    return p.parents[0]

# collect paths to things we need
gccbin = findgccbin()
gcclib64 = gccbin.parents[0].joinpath('lib64')
if not gcclib64.is_dir():
    gcclib64 = gccbin.parents[0].joinpath('lib')

# define include dirs
GCCLIB64 = gcclib64.as_posix()

ldirs=['src/', '/usr/local/lib', GCCLIB64]
link_list=['-fopenmp']

EXTENSIONS = Extension(
    name="qihmc.qihmc_cc",
    sources=SOURCES,
    library_dirs = ldirs,
    include_dirs=[SRC, get_include()],
    language="c++14",
    extra_compile_args=["-O0", "-march=native", "-ffast-math", "-fopenmp"],
    extra_link_args = link_list,
)

with open(path.join(ROOT, "README.md"), encoding="utf-8") as inp:
    LONG_DESCRIPTION = inp.read()

setup(
    name="qihmc",
    version=__version__,
    description=None,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url=None,
    author=__author__,
    author_email=None,
    keywords=[],
    packages=find_packages(exclude=["docs", "tests"]),
    package_dir={'': 'src'},
    ext_modules=cythonize([EXTENSIONS]),
)
