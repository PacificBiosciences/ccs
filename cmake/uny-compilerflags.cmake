
include(CheckCXXCompilerFlag)

# shared CXX flags for all source code & tests
set(UNY_FLAGS "-std=c++11 -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable")

# gperftools support
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(UNY_LINKER_FLAGS "${UNY_LINKER_FLAGS} -Wl,-no_pie")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

# NOTE: quash clang warnings w/ Boost
check_cxx_compiler_flag("-Wno-unused-local-typedefs" HAS_NO_UNUSED_LOCAL_TYPEDEFS)
if(HAS_NO_UNUSED_LOCAL_TYPEDEFS)
    set(UNY_FLAGS "${UNY_FLAGS} -Wno-unused-local-typedefs")
endif()

# Coverage settings
set(UNY_DEBUG_FLAGS "${UNY_FLAGS} -fprofile-arcs -ftest-coverage")

# Extra testing that will lead to longer compilation times!
if (SANITIZE)
    # AddressSanitizer is a fast memory error detector
    set(UNY_DEBUG_FLAGS "${UNY_DEBUG_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")

    # Clang Thread Safety Analysis is a C++ language extension which warns about 
    # potential race conditions in code.
    set(UNY_DEBUG_FLAGS "${UNY_DEBUG_FLAGS} -Wthread-safety")

    # ThreadSanitizer is a tool that detects data races
    set(UNY_DEBUG_FLAGS "${UNY_DEBUG_FLAGS} -fsanitize=thread")

    # MemorySanitizer is a detector of uninitialized reads.
    set(UNY_DEBUG_FLAGS "${UNY_DEBUG_FLAGS} -fsanitize=memory")

    # UndefinedBehaviorSanitizer is a fast undefined behavior detector.
    set(UNY_DEBUG_FLAGS "${UNY_DEBUG_FLAGS} -fsanitize=undefined")
    MESSAGE(${UNY_DEBUG_FLAGS})
endif()