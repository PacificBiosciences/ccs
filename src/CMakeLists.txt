include_directories (${Boost_INCLUDE_DIRS}
                     )
add_executable (ChimeraLabeler ChimeraLabeler.cpp)
target_link_libraries (ChimeraLabeler
                       ${Boost_FILESYSTEM_LIBRARY}
                       ${Boost_SYSTEM_LIBRARY}
                       ${ZLIB_LIBRARIES}
                       )
