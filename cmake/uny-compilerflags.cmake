
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# shared CXX flags for all source code & tests
set(UNY_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-unused-variable")

# gperftools support
if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND APPLE)
    set(UNY_LINKER_FLAGS "${UNY_LINKER_FLAGS} -Wl,-no_pie")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug" AND APPLE)

# static linking
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(UNY_LINKER_FLAGS "${UNY_LINKER_FLAGS} -static-libgcc -static-libstdc++")
ENDIF()

# NOTE: quash clang warnings w/ Boost
check_cxx_compiler_flag("-Wno-unused-local-typedefs" HAS_NO_UNUSED_LOCAL_TYPEDEFS)
if(HAS_NO_UNUSED_LOCAL_TYPEDEFS)
    set(UNY_FLAGS "${UNY_FLAGS} -Wno-unused-local-typedefs")
endif()

# Cannot use this until pbbam complies
# if (CMAKE_COMPILER_IS_GNUCXX)
#     set(UNY_FLAGS "${UNY_FLAGS} -Werror=suggest-override")
# endif()

# Coverage settings
if (UNY_inc_coverage)
    set(UNY_COV_FLAGS "${UNY_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

# Extra testing that will lead to longer compilation times!
if (SANITIZE)
    # AddressSanitizer is a fast memory error detector
    set(UNY_SANITY_FLAGS "${UNY_SANITY_FLAGS} -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")

    # Clang Thread Safety Analysis is a C++ language extension which warns about
    # potential race conditions in code.
    set(UNY_SANITY_FLAGS "${UNY_SANITY_FLAGS} -Wthread-safety")

    # ThreadSanitizer is a tool that detects data races
    set(UNY_SANITY_FLAGS "${UNY_SANITY_FLAGS} -fsanitize=thread")

    # MemorySanitizer is a detector of uninitialized reads.
    set(UNY_SANITY_FLAGS "${UNY_SANITY_FLAGS} -fsanitize=memory")

    # UndefinedBehaviorSanitizer is a fast undefined behavior detector.
    set(UNY_SANITY_FLAGS "${UNY_SANITY_FLAGS} -fsanitize=undefined")
endif()

# shared CXX flags for src & tests
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UNY_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UNY_COV_FLAGS} ${UNY_SANITY_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${UNY_LINKER_FLAGS}")
