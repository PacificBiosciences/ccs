// Author: Lance Hepler

#pragma once

#include <cstdio>
#include <cstdlib>

#include <sys/ioctl.h>

namespace PacBio {
namespace Util {

inline void SetColumns()
{
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    static char env[32];
    unsigned short col_size = w.ws_col > 0 ? w.ws_col : 80;
    int retval = snprintf(env, sizeof(env), "COLUMNS=%u", col_size);
    if (retval >= 0 && static_cast<size_t>(retval) < sizeof(env)) {
        putenv(env);
    }
}

}  // namespace Util
}  // namespace PacBio
