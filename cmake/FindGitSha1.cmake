
if(__find_git_sha1)
    return()
endif()
set(__find_git_sha1 YES)

function(find_git_sha1 _GIT_SHA1)
    find_package(Git QUIET REQUIRED)
    execute_process(COMMAND
        "${GIT_EXECUTABLE}" "describe" "--always" "--dirty=*"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        RESULT_VARIABLE res
        OUTPUT_VARIABLE out
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (NOT res EQUAL 0)
        message(FATAL_ERROR "Could not determine git sha1 via `git describe --always --dirty=*`")
    endif()
    set(${_GIT_SHA1} "${out}" PARENT_SCOPE)
endfunction()
