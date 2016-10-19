
from __future__ import print_function

import os
import os.path
import re
import sys

from copy import copy
from distutils.command.build_ext import build_ext
from distutils.util import strtobool
from setuptools import setup, Extension
from shutil import copy2, rmtree
from subprocess import Popen

def ParseVersion():
    thisDir = os.path.dirname(os.path.realpath(__file__))
    cmakeLists = os.path.join(thisDir, "CMakeLists.txt")
    regexp = re.compile(r'project\([^ ]+ VERSION (\d+\.\d+\.\d+) [^\)]+\)')
    with open(cmakeLists) as handle:
        for line in handle:
            m = regexp.search(line)
            if m:
                return m.group(1)
    raise RuntimeError("unable to find version string!")

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
        value = self.env.get(key, default)
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

    def __call__(self, sourceDir, binaryDir, targets=[]):
        configure = copy(self.configure)
        if self.generator:
            configure.append("-G{0}".format(self.generator))
        for k, v in self.definitions.iteritems():
            configure.append("-D{0}={1}".format(k, v))
        configure.append(sourceDir)
        print("Configuring with command `{0}`".format(" ".join(configure)), file=sys.stderr)
        retcode = Popen(configure, cwd=binaryDir, env=self.env).wait()
        if retcode != 0:
            raise RuntimeError("failed to configure the project")
        build = self.build + targets
        print("Building with command `{0}`".format(" ".join(build)), file=sys.stderr)
        retcode = Popen(build, cwd=binaryDir).wait()
        if retcode != 0:
            raise RuntimeError("failed to build the project")

class MyBuildExt(build_ext):
    def build_extension(self, ext):
        destDir = os.path.dirname(self.get_ext_fullpath(ext.name))
        thisDir = os.path.dirname(os.path.realpath(__file__))
        try:
            os.makedirs(self.build_temp)
            os.makedirs(destDir)
        except OSError:
            pass
        env = os.environ.copy()
        cmake = CMake(env)
        cmake.add_definition_from_env("Boost_INCLUDE_DIRS")
        cmake.add_definition_from_env("PYTHON_INCLUDE_DIRS")
        cmake.add_definition_from_env("pbcopper_INCLUDE_DIRS")
        cmake.add_definition_from_env("pbcopper_LIBRARIES")
        cmake.add_definition_from_env("GIT_EXECUTABLE")
        cmake.add_definition_from_env("SWIG_COMMAND")
        cmake.add_definition_from_env("UNY_use_ccache")
        cmake.add_definition_from_env("CMAKE_BUILD_TYPE", "RelWithDebInfo")
        cmake.add_definition("PYTHON_SWIG", "1")
        cmake.add_definition("UNY_build_tests", "0")
        cmake.add_definition("UNY_build_bin", "0")
        cmake.add_definition("PYTHON_EXECUTABLE", sys.executable)
        cmake.set_generator("Ninja" if which("ninja", env) else "Default")
        targets = ["_ConsensusCore2"]
        try:
            cmake(thisDir, self.build_temp, targets)
        except RuntimeError:
            # if this happens we're probably under pip and need to regenerate
            rmtree(self.build_temp)
            os.makedirs(self.build_temp)
            cmake(thisDir, self.build_temp, targets)
        for fname in ["_ConsensusCore2.so", "ConsensusCore2.py"]:
            copy2(os.path.join(self.build_temp, "swig", "lib", fname), destDir)

setup(
    name="ConsensusCore2",
    version=ParseVersion(),
    author="PacificBiosciences",
    author_email="devnet@pacb.com",
    url="http://www.github.com/PacificBiosciences/ConsensusCore2",
    description="A library for generating consensus sequences for PacBio data",
    license="BSD",
    ext_modules=[Extension("_ConsensusCore2", [])],
    package_dir={"": "swig"},
    py_modules=["ConsensusCore2"],
    install_requires=["numpy >= 1.6.0"],
    setup_requires=["numpy >= 1.6.0"],
    cmdclass={"build_ext": MyBuildExt}
    )
