CONFIG_MAKE  = make/Config.mk
ifneq ("$(wildcard $(CONFIG_MAKE))","")
include $(CONFIG_MAKE)
else
$(error Run ./configure first!)
endif

include make/Defs.mk

ifeq (,$(wildcard $(PYTHON_INCLUDE)/Python.h))
$(error python include not configured properly, cannot build python module!)
endif

ifeq (,$(wildcard $(NUMPY_INCLUDE)/numpy/arrayobject.h))
$(error numpy include not configured properly, cannot build python module!)
endif

SWIG_CMD   := SWIG_LIB=$(SWIG_LIB) $(SWIG) -Wall -c++ -python -builtin
GEN_CXX    := $(PYTHON_BUILD_DIR)/ConsensusCore_wrap.cxx
PYTHON_DLL := $(PYTHON_BUILD_DIR)/_ConsensusCore.so

all: $(PYTHON_DLL)

$(PYTHON_DLL): $(SWIG_INTERFACES) $(CXX_LIB)
	-mkdir -p $(PYTHON_BUILD_DIR)
	$(SWIG_CMD) $(INCLUDES) -module ConsensusCore -o $(GEN_CXX) $(SWIG_INTERFACE)
	$(CXX) $(SHLIB_FLAGS) $(INCLUDES) -I $(PYTHON_INCLUDE) -I $(NUMPY_INCLUDE) $(GEN_CXX) $(CXX_LIB) -o $(PYTHON_DLL)

test-python: $(PYTHON_DLL)
	@PYTHONPATH=$(PYTHON_BUILD_DIR) python src/Demos/Demo.py && echo "Python build is OK!"

.PHONY: all test-python $(PYTHON_DLL)
