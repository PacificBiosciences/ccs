#
# Basic makefile. Use this for building the C++ library and the SWIG
# bindings.
#
# Notes:
#  - can compile with clang++ using:
#      % make GXX=clang++
#  - 32-bit build:
#      % make MACHINE=-m32
#  - debug build:
#      % make DEBUG=1
#
include make/Config.mk
include make/Defs.mk

all: lib

lib:
	${MAKE} -f make/Cpp.mk

#
# SWIG targets
#
python: lib
	$(MAKE) -f make/Python.mk

echo-python-build-directory:
	@echo $(PYTHON_BUILD_DIR)

csharp: lib
	$(MAKE) -f make/CSharp.mk

before-xbuild: csharp
	-mkdir -p bin/Debug
	-mkdir -p bin/Release
	-cp $(BUILD_ROOT)/CSharp/libConsensusCore.so bin/Debug/
	-cp $(BUILD_ROOT)/CSharp/libConsensusCore.so bin/Release/

#
# Clean targets
#
clean-python:
	-rm -rf $(PYTHON_BUILD_DIR)
	-rm -rf ConsensusCore.egg-info
	-rm -rf dist

clean-csharp:
	-rm -rf $(CSHARP_BUILD_DIR)
	-rm -rf bin obj

clean-cxx:
	-rm -rf $(BUILD_ROOT)/C++

clean: clean-cxx clean-python clean-csharp
	-rm -rf $(BUILD_ROOT)

#
# Test targets
#
test-python:
	@make -f make/Python.mk test-python

test-csharp:
	@make -f make/CSharp.mk test-csharp

test: lib
	@make -f make/Tests.mk

check: test
tests: test

arrowcheck: lib
	$(CXX) $(CXX_FLAGS) $(INCLUDES) src/Demos/MatrixTester.cpp $(OBJDIR)/libConsensusCore.a -o $(OBJDIR)/arrowcheck
	$(OBJDIR)/arrowcheck

#
# Lint targets
#

lint:
	-find src -name "*.[ch]pp" | xargs ./tools/cpplint.py --verbose=0 --counting=toplevel

pre-commit-hook:
#       for speed, apply cpplint only to changed files.
	git diff --cached --name-only --diff-filter=ACM | \
         grep -e  '.*.[ch]pp$$' | xargs tools/cpplint.py --verbose=3


#
# Targets used by PBI internal build
#
pip-uninstall: $(shell which pip > /dev/null)
	@pip freeze|grep 'ConsensusCore=='>/dev/null \
      && pip uninstall -y ConsensusCore \
      || true

pip-install: $(shell which pip > /dev/null)
	@pip install --no-index \
          --install-option="--swig=$(SWIG)" \
          --install-option="--boost=$(BOOST)" \
          ./

.PHONY: all lib clean-cxx clean test tests check python clean-python \
	csharp clean-csharp echo-python-build-directory \
	test-python test-csharp pip-uninstall pip-install \
	lint pre-commit-hook 
