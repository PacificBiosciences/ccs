
from __future__ import print_function

import os
import os.path
import sys

from distutils import sysconfig
from distutils.command.clean import clean as dclean
from distutils.errors import CompileError
from setuptools import setup
from setuptools.dist import Distribution
from shutil import rmtree
from subprocess import Popen
from tempfile import mkdtemp

thisDir = os.path.dirname(os.path.realpath(__file__))
swigLib = os.path.join(thisDir, "swig", "lib")
modName = "_ConsensusCore2.so"

class MyClean(dclean):
    def run(self):
        rmtree(swigLib, ignore_errors=True)
        dclean.run(self)

class MyDist(Distribution):
    def has_ext_modules(self):
        return True

if not os.path.exists(os.path.join(swigLib, modName)):
    buildDir = mkdtemp()
    try:
        env = os.environ.copy()
        cmake = env.get("CMAKE_COMMAND", "cmake")
        boost = env.get("Boost_INCLUDE_DIRS", None)
        pyinc = env.get("PYTHON_INCLUDE_DIRS", None)
        rpath = env.get("CMAKE_SKIP_RPATH", None)
        swig = env.get("SWIG_COMMAND", None)
        cmds = [cmake,
                "-DCMAKE_BUILD_TYPE=RelWithDebInfo",
                "-DPYTHON_SWIG=1",
                "-DPYTHON_EXECUTABLE={0}".format(sys.executable)]
        if boost is not None:
            cmds.append("-DBoost_INCLUDE_DIRS={0}".format(boost))
        if pyinc is not None:
            cmds.append("-DPYTHON_INCLUDE_DIRS={0}".format(pyinc))
        if rpath is not None:
            cmds.append("-DCMAKE_SKIP_RPATH=TRUE")
        if swig is not None:
            cmds.append("-DSWIG_COMMAND={0}".format(swig))
        cmds.append(thisDir)
        print("Running command {0}".format(" ".join(cmds)), file=sys.stderr)
        retcode = Popen(cmds, cwd=buildDir, env=env).wait()
        if (retcode != 0):
            raise RuntimeError("failed to configure the ConsensusCore2 with CMake!")
        verbose = env.get("VERBOSE", None)
        cmds = ["make"]
        if verbose is not None:
            cmds.append("VERBOSE=1")
        print("Running command {0}".format(" ".join(cmds)), file=sys.stderr)
        retcode = Popen(cmds, cwd=buildDir).wait()
        if (retcode != 0):
            raise RuntimeError("failed to compile or link ConsensusCore2!")
    finally:
        rmtree(buildDir)

setup(
    name="ConsensusCore2",
    version="0.11.0",
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
