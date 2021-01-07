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

    # determine if we have clang or gnu gcc or raise an error
    have_gnu = any(['gcc version' in p for p in ps])
    have_clang = any(['clang' in p for p in ps])

    if not have_gnu and not have_clang:
        raise SystemError("Unsupported compiler.  We only support gnu gcc and clang now")

    p = pathlib.Path(gcc)
    # set some associated env vars for make
    environ["CC"] = p.as_posix()
    # assume that CXX is in the same directory as CC and has the same extension
    # for example, c++-10 -> g++-10, gcc-10 -> g++-10
    gcc_path = str(pathlib.Path(environ["CC"]).parent.absolute())
    gxx      = environ["CC"].split('/')[-1].replace('gcc','g++').replace('c++','g++')
    environ["CXX"] = gcc_path +'/'+ gxx
    if not path.exists(environ["CXX"]):
        raise SystemError("Missing g++")
    else:
        print("Found g++ at ",environ["CXX"])
    # Fix failure on unbuntu
    # https://stackoverflow.com/questions/45308426/distutils-partially-ignoring-cc-environment-variable-when-building-python-extens
    if platform.system() != 'Darwin':
        environ["LDSHARED"] = environ["CXX"] + " -std=c++14 -pthread -shared"
    return p.parents[0], have_clang

# collect paths to things we need
gccbin, have_clang = findgccbin()
gcclib64 = gccbin.parents[0].joinpath('lib64')
if not gcclib64.is_dir():
    gcclib64 = gccbin.parents[0].joinpath('lib')

# define include dirs
GCCLIB64 = gcclib64.as_posix()
ldirs=['src/', '/usr/local/lib', GCCLIB64]
i_dirs = [SRC, get_include()]
if have_clang:
    link_list    = ['-Xpreprocessor', '-fopenmp', '-stdlib=libc++']
    compile_args = ["-O0", "-march=native", "-ffast-math", "-Xpreprocessor", "-fopenmp", "-stdlib=libc++", "-std=c++14"]
    if 'Xcode' in str(gccbin):
        #include_dir = '/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include'
        raise SystemError("We do not yet know how to point to the actual clang compiler in Xcode.app\nTry pointing to either\n/usr/bin/c++\n/usr/bin/gcc\nif you want to use the apple clang")
else:
    link_list    = ['-fopenmp']
    compile_args = ["-O0", "-march=native", "-ffast-math", "-fopenmp"]

EXTENSIONS = Extension(
    name="qihmc.qihmc_cc",
    sources=SOURCES,
    library_dirs = ldirs,
    include_dirs = i_dirs,
    language="c++14",
    extra_compile_args=compile_args,
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
