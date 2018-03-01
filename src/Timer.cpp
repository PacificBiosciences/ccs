// Authors: Lance Hepler, Armin Töpfer

#include <chrono>
#include <sstream>

#include <pacbio/util/Timer.h>

using std::chrono::steady_clock;
using std::chrono::milliseconds;

namespace PacBio {
namespace Util {

Timer::Timer() { Restart(); }

float Timer::ElapsedMilliseconds() const
{
    auto tock = steady_clock::now();
    return std::chrono::duration_cast<milliseconds>(tock - tick).count();
}

float Timer::ElapsedSeconds() const { return ElapsedMilliseconds() / 1000; }

std::string Timer::ElapsedTime() const
{
    auto tock = steady_clock::now();
    auto t = std::chrono::duration_cast<std::chrono::nanoseconds>(tock - tick).count();

    auto d = t / 10000 / 1000 / 1000 / 60 / 60 / 24;
    auto h = (t / 1000 / 1000 / 1000 / 60 / 60) % 24;
    auto m = (t / 1000 / 1000 / 1000 / 60) % 60;
    auto s = (t / 1000 / 1000 / 1000) % 60;
    auto ms = (t / 1000 / 1000) % 1000;
    auto us = (t / 1000) % 1000;
    auto ns = t % 1000;
    std::ostringstream ss;
    if (d > 0) ss << d << "d ";
    if (h > 0) ss << h << "h ";
    if (m > 0) ss << m << "m ";
    if (s > 0) ss << s << "s ";
    if (ms > 0) ss << ms << "ms ";
    if (us > 0) ss << us << "us ";
    if (ns > 0) ss << ns << "ns ";
    return ss.str();
}

void Timer::Restart() { tick = steady_clock::now(); }

}  // namespace Util
}  // namespace PacBio
