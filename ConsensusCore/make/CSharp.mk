include make/Config.mk
include make/Defs.mk

SWIG_CMD        = SWIG_LIB=$(SWIG_LIB) $(SWIG) -Wall -c++ -csharp

CSHARP_BUILD_DIR = build/CSharp
GEN_CXX          = $(CSHARP_BUILD_DIR)/ConsensusCore_wrap.cxx
PINVOKE_LIB      = $(CSHARP_BUILD_DIR)/libConsensusCore.so

$(PINVOKE_LIB): $(SWIG_INTERFACES) $(CXX_LIB)
	mkdir -p $(CSHARP_BUILD_DIR)
	$(SWIG_CMD) $(INCLUDES) -module ConsensusCore -namespace ConsensusCore \
	        -o $(GEN_CXX) -outdir $(CSHARP_BUILD_DIR) $(SWIG_INTERFACE)
	cat $(CSHARP_BUILD_DIR)/*.cs > $(CSHARP_BUILD_DIR)/.ConsensusCore.cs
	rm $(CSHARP_BUILD_DIR)/*.cs
	mv $(CSHARP_BUILD_DIR)/.ConsensusCore.cs $(CSHARP_BUILD_DIR)/ConsensusCore.cs
	$(CXX) $(SHLIB_FLAGS) $(INCLUDES) $(GEN_CXX) $(CXX_LIB) -o $(PINVOKE_LIB)

test-csharp: $(PINVOKE_LIB)
	@(cd $(CSHARP_BUILD_DIR) && \
	 mcs -target:library ConsensusCore.cs && \
	 mcs $(PROJECT_ROOT)/src/Demos/Demo.cs -reference:ConsensusCore.dll -out:Demo.exe && \
	 mono Demo.exe && \
	 echo "CSharp build OK!")

.PHONY: all clean test-csharp
