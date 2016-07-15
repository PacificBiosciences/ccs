
from __future__ import print_function

import os
import os.path
import sys

from copy import copy
from distutils import sysconfig
from distutils.command.clean import clean as dclean
from distutils.errors import CompileError
from distutils.util import strtobool
from setuptools import setup
from setuptools.dist import Distribution
from shutil import rmtree
from subprocess import Popen
from tempfile import mkdtemp

thisDir = os.path.dirname(os.path.realpath(__file__))
swigLib = os.path.join(thisDir, "swig", "lib")
modName = "_ConsensusCore2.so"

def which(program, env=os.environ):
    def isExe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if isExe(program):
            return program
    else:
        for path in env.get("PATH", "").split(os.pathsep):
            path = path.strip("\"")
            exe = os.path.join(path, program)
            if isExe(exe):
                return exe
    return None


class MyClean(dclean):
    def run(self):
        rmtree(swigLib, ignore_errors=True)
        dclean.run(self)


class MyDist(Distribution):
    def has_ext_modules(self):
        return True


class CMake(object):
    def __init__(self, env=os.environ, verbose=False):
        cmake = env.get("CMAKE_COMMAND", "cmake") or which("cmake", env)
        if cmake is None:
            raise RuntimeError("cannot find `cmake` command, "
                    "please populate CMAKE_COMMAND environment variable")
        self.build = ["make"]
        self.env = env
        self.configure = [cmake]
        self.definitions = dict()
        self.generator = None

        rpath = self.env.get("CMAKE_SKIP_RPATH", "")
        try:
            rpath = bool(strtobool(rpath))
        except ValueError:
            rpath = False
        if rpath:
            self.add_definition("CMAKE_SKIP_RPATH", "TRUE")

        verbose = self.env.get("VERBOSE", "{0}".format(verbose))
        try:
            verbose = bool(strtobool(verbose))
        except ValueError:
            verbose = False
        self.verbose = verbose

    def add_definition(self, key, value):
        self.definitions[key] = value

    def add_definition_from_env(self, key, default=None):
        value = env.get(key, default)
        if value:
            self.add_definition(key, value)

    def set_build_type(self, buildType):
        if buildType not in set(["Release", "Debug", "RelWithDebInfo"]):
            raise ValueError("CMAKE_BUILD_TYPE must be in (Release, Debug, RelWithDebInfo)")
        self.add_definition("CMAKE_BUILD_TYPE", buildType)

    def set_generator(self, gen):
        if gen not in set(["Default", "Ninja"]):
            raise ValueError("valid generators must be in (default, ninja)")
        if gen == "Default":
            gen = None
            self.build = ["make"]
            if self.verbose:
                self.build.append("VERBOSE=1")
        elif gen == "Ninja":
            self.ninja = self.env.get("NINJA_COMMAND", None) or which("ninja", self.env)
            if not self.ninja:
                raise RuntimeError("cannot find `ninja` command, "
                        "please populate NINJA_COMMAND environment variable")
            self.build = [self.ninja]
            if self.verbose:
                self.build.append("-v")
        self.generator = gen

    def __call__(self, directory):
        configure = copy(self.configure)
        if self.generator:
            configure.append("-G{0}".format(self.generator))
        for k, v in self.definitions.iteritems():
            configure.append("-D{0}={1}".format(k, v))
        configure.append(directory)
        print("Configuring with command `{0}`".format(" ".join(configure)), file=sys.stderr)
        retcode = Popen(configure, cwd=buildDir, env=env).wait()
        if retcode != 0:
            raise RuntimeError("failed to configure the project")
        print("Building with command `{0}`".format(" ".join(self.build)), file=sys.stderr)
        retcode = Popen(self.build, cwd=buildDir).wait()
        if retcode != 0:
            raise RuntimeError("failed to build the project")


# build the library if we haven't already
if not os.path.exists(os.path.join(swigLib, modName)):
    buildDir = mkdtemp()
    try:
        env = os.environ.copy()
        cmake = CMake(env)
        cmake.add_definition_from_env("Boost_INCLUDE_DIRS")
        cmake.add_definition_from_env("PYTHON_INCLUDE_DIRS")
        cmake.add_definition_from_env("SWIG_COMMAND")
        cmake.add_definition("PYTHON_SWIG", "1")
        cmake.add_definition("PYTHON_EXECUTABLE", sys.executable)
        cmake.set_build_type("RelWithDebInfo")
        cmake.set_generator("Ninja" if which("ninja", env) else "Default")
        cmake(thisDir)
    finally:
        rmtree(buildDir)

setup(
    name="ConsensusCore2",
    version="0.13.0",
    author="PacificBiosciences",
    author_email="devnet@pacb.com",
    url="http://www.github.com/PacificBiosciences/ConsensusCore2",
    description="A library for generating consensus sequences for PacBio data",
    license="BSD",
    packages=["ConsensusCore2"],
    package_dir={"ConsensusCore2": swigLib},
    package_data={"ConsensusCore2": [modName]},
    install_requires=["numpy >= 1.6.0"],
    setup_requires=["numpy >= 1.6.0"],
    distclass=MyDist,
    cmdclass={"clean": MyClean}
    )
