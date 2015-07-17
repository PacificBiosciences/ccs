THIS_DIRECTORY := $(dir $(lastword $(MAKEFILE_LIST)))
PROJECT_ROOT   := $(abspath $(THIS_DIRECTORY)/..)
BUILD_ROOT     ?= $(PROJECT_ROOT)/build
GMOCK_ROOT     ?= $(PROJECT_ROOT)/external/gmock-1.7.0

# Without this, rm -rf *.(o,so) fails
SHELL=/bin/bash

INCLUDES        := -I$(PROJECT_ROOT)/include -I$(PROJECT_ROOT)/external

ifdef COVERAGE
OBJDIR          := $(BUILD_ROOT)/C++/Coverage
else
OBJDIR          := $(BUILD_ROOT)/C++
endif

VPATH           := src/C++/                     \
                  :src/C++/Align                \
                  :src/C++/Arrow                \
                  :src/C++/Arrow/Matrix         \
                  :src/C++/Arrow/detail         \
                  :src/C++/Edna                 \
                  :src/C++/Matrix/              \
                  :src/C++/Quiver/              \
                  :src/C++/Quiver/detail        \
                  :src/C++/Poa/                 \
                  :src/C++/Statistics           \
                  :src/C++/Logging

CXX_LIB         := $(abspath $(OBJDIR)/libConsensusCore.a)
CXX_SRCS        := $(patsubst src/C++/%,%,$(shell find src/C++ -name "*.cpp" | grep -v '\#'))

CXX_OPT_FLAGS_DEBUG   := -O0 -g
CXX_OPT_FLAGS_RELEASE := -O3 -DNDEBUG -g

ifeq ($(DEBUG),)
        CXX_OPT_FLAGS = $(CXX_OPT_FLAGS_RELEASE)
else
        CXX_OPT_FLAGS = $(CXX_OPT_FLAGS_DEBUG)
endif

# Detect mac/linux
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
        GXX ?= clang++
else
        GXX ?= g++
endif

ifeq ($(GXX),clang++)
    CXX_FLAGS           = $(GXX_FLAGS) $(CXX_OPT_FLAGS) -msse3 -fPIC -Qunused-arguments -fno-omit-frame-pointer
    CXX_STRICT_FLAGS    = $(GXX_FLAGS) $(CXX_FLAGS) -pedantic -Wall
else
    CXX_FLAGS           = $(CXX_OPT_FLAGS) $(CXX_EXTRA_ARGS) -msse3 -fPIC -fno-omit-frame-pointer
    CXX_STRICT_FLAGS    = $(CXX_FLAGS) -pedantic -Wall
endif

CXX             = $(CCACHE) $(GXX) $(MACHINE) -std=$(CPP_ABI) $(CXX_FLAGS) $(INCLUDES) -isystem $(BOOST)
CXX_STRICT      = $(CCACHE) $(GXX) $(MACHINE) -std=$(CPP_ABI) $(CXX_STRICT_FLAGS) $(INCLUDES) -isystem $(BOOST)

ifeq ($(UNAME), Darwin)
    SHLIB_FLAGS = -shared -undefined dynamic_lookup
else
    SHLIB_FLAGS = -pthread -shared -Wl,-O1
endif

GMOCK_LIBSRC := $(GMOCK_ROOT)/gmock-gtest-all.cc
GMOCK_MAIN   := $(GMOCK_ROOT)/gmock_main.cc

PYTHON_BUILD_DIR := $(BUILD_ROOT)/Python
CSHARP_BUILD_DIR := $(BUILD_ROOT)/CSharp

SWIG_INTERFACE  = src/SWIG/ConsensusCore.i
SWIG_INTERFACES = $(shell find src/SWIG/ -name "*.i" | grep -v '\#')
