#!/usr/bin/env python
from __future__ import print_function

from distutils.command.build import build as _build
from distutils.errors import CompileError
from distutils import sysconfig
from setuptools import setup

import os, re, sys
from glob import glob
from os.path import join, dirname as up

def die(msg):
    print(msg)
    sys.exit(1)

def str_version(tpl):
    return ".".join(map(str, tpl))

def configure():
    """
    Look up the dependencies
     - Python.h
     - numpy headers
     - boost headers (from configuration flag)

    meanwhile, pass through any config-looking arguments (ex: --boost=...)
    """
    # Python.h
    python_inc = sysconfig.get_python_inc()

    # Libpython
    python_libdir = join(sysconfig.get_python_lib(standard_lib=True), "config")
    python_lib = "python" + sysconfig.get_python_version()

    if '--version' in sys.argv:
        return "--version"
    # numpy headers
    try:
        import numpy
        numpy_inc = numpy.get_include()
    except ImportError:
        die("Requires numpy >= 1.6.0")

    ccArgs = re.compile(r"^--(?:boost|swig|swig-lib)=?")
    ccOpts = re.compile(r"^--(?:debug|c\+\+11|pbi|modules)$")

    configArgs = []
    for arg in sys.argv[:]:
        if ccArgs.match(arg) or ccOpts.match(arg):
            configArgs.append(arg)
            sys.argv.remove(arg)

    configArgs.append("--python-include=%s " % python_inc)
    configArgs.append("--numpy-include=%s "  % numpy_inc)

    return " ".join(configArgs)

# This always has to be run, to remove extra arguments that will
# confuse setuptools otherwise
configuration = configure()


class build(_build):
    """
    Build the native code using the Makefile, sidestepping
    most of the horrors of Python packaging.
    """
    def run(self):
        error = os.system("./configure " + configuration)
        if error:
            raise CompileError("Failed to configure ConsensusCore build")
        error = os.system("make python")
        if error:
            raise CompileError("Failed to compile or link ConsensusCore C++ code")

    @staticmethod
    def pythonBuildDirectory():
        """
        Returns the directory where the build stashed the generated
        Python module.
        """
        #return os.popen("make --no-print-directory echo-python-build-directory").read().strip()
        return "build/Python"

# HACK: "setup.py install" will fail if build hasn't been called yet,
# because setuptools attempts to scan the package directory before
# invoking the build---which will fail because the directory is
# created by the build.  I consider this a bug in setuptools, and this
# is the least involved workaround I have found: inject "build" into
# the command line.
if "install" in sys.argv and not "build" in sys.argv:
    installPos = sys.argv.index("install")
    sys.argv.insert(installPos, "build")

setup(name="ConsensusCore",
      version="1.0.1",
      author="Pacific Biosciences",
      author_email="devnet@pacificbiosciences.com",
      url="http://www.github.com/PacificBiosciences/ConsensusCore",
      description= \
          """A library for genomic consensus and variant calling""",
      license=open("LICENSES").read(),
      py_modules=["ConsensusCore"],
      packages = [""],
      package_dir={"": build.pythonBuildDirectory() },

      # Smuggle the native library in as a data file
      package_data={"": ["_ConsensusCore.so"]},

      cmdclass={"build": build})
