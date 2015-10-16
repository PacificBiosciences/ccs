
from __future__ import print_function

import os
import os.path
import sys

from distutils import sysconfig
from distutils.command.build import build as dbuild
from distutils.command.build_py import build_py as dbuild_py
from distutils.command.clean import clean as dclean
from distutils.command.install import install as dinstall
from distutils.errors import CompileError
from setuptools import setup
from setuptools.dist import Distribution
from shutil import rmtree
from subprocess import Popen
from tempfile import mkdtemp

thisDir = os.path.dirname(os.path.realpath(__file__))
swigLib = os.path.join(thisDir, "swig", "lib")

def cleanAll(self):
    cleaner = dclean(self.distribution)
    cleaner.all = True
    cleaner.finalize_options()
    cleaner.run()
    rmtree(swigLib, ignore_errors=True)

class MyBuild(dbuild_py):
    def run(self):
        buildDir = mkdtemp()
        try:
            cleanAll(self)
            env = os.environ.copy()
            cmake = env.get("CMAKE_COMMAND", "cmake")
            boost = env.get("Boost_INCLUDE_DIRS", None)
            swig = env.get("SWIG_COMMAND", None)
            cmds = [cmake,
                    "-DCMAKE_BUILD_TYPE=RelWithDebInfo",
                    "-DPYTHON_SWIG=1",
                    "-DPYTHON_EXECUTABLE={0}".format(sys.executable)]
            if boost is not None:
                cmds.append("-DBoost_INCLUDE_DIRS={0}".format(boost))
            if swig is not None:
                cmds.append("-DSWIG_COMMAND={0}".format(swig))
            cmds.append(thisDir)
            print("running command: {0}".format(" ".join(cmds)), file=sys.stderr)
            retcode = Popen(cmds, cwd=buildDir, env=env).wait()
            if (retcode != 0):
                raise CompileError("failed to configure the pbconsensus with CMake!")
            print("running command: make", file=sys.stderr)
            retcode = Popen(["make"], cwd=buildDir).wait()
            if (retcode != 0):
                raise CompileError("failed to compile or link pbconsensus!")
        finally:
            rmtree(buildDir)
            try: super(MyBuild, self).run()
            except TypeError: dbuild_py.run(self)

class MyInstall(dinstall):
    def run(self):
        try: super(MyInstall, self).run()
        except TypeError: dinstall.run(self)
        finally: cleanAll(self)

class MyDist(Distribution):
    def has_ext_modules(self):
        return True

setup(
    name="pbconsensus",
    version="0.9.0",
    author="PacificBiosciences",
    author_email="devnet@pacb.com",
    url="http://www.github.com/PacificBiosciences/pbconsensus",
    description="A library for generating consensus sequences for PacBio data",
    license="BSD",
    packages=["pbconsensus"],
    package_dir={"pbconsensus": swigLib},
    package_data={"pbconsensus": ["_pbconsensus.so"]},
    install_requires=["numpy >= 1.6.0"],
    setup_requires=["numpy >= 1.6.0"],
    distclass=MyDist,
    cmdclass={"build_py": MyBuild}
    )
