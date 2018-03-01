// Author: Lance Hepler

#pragma once

#include <chrono>
#include <string>

namespace PacBio {
namespace Util {

class Timer
{
public:
    Timer();

    float ElapsedMilliseconds() const;
    float ElapsedSeconds() const;
    std::string ElapsedTime() const;
    void Restart();

private:
    std::chrono::time_point<std::chrono::steady_clock> tick;
};

}  // namespace Util
}  // namespace PacBio
