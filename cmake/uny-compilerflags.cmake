
include(CheckCXXCompilerFlag)

# shared CXX flags for all source code & tests
if (MSVC)
    set(UNY_FLAGS "/Wall")
else()
    set(UNY_FLAGS "-std=c++11 -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable")
endif()

# NOTE: quash clang warnings w/ Boost
check_cxx_compiler_flag("-Wno-unused-local-typedefs" HAS_NO_UNUSED_LOCAL_TYPEDEFS)
if(HAS_NO_UNUSED_LOCAL_TYPEDEFS)
    set(UNY_FLAGS "${UNY_FLAGS} -Wno-unused-local-typedefs")
endif()