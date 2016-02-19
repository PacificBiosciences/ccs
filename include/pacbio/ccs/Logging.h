// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Lance Hepler

#pragma once

#include <chrono>
#include <condition_variable>
#include <csignal>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

namespace PacBio {
namespace CCS {
namespace Logging {

enum LogLevel : unsigned char
{
    TRACE = 0,
    DEBUG = 10,
    INFO = 20,
    NOTICE = 30,
    WARN = 40,
    ERROR = 50,
    CRITICAL = 60,
    FATAL = 70
};

inline LogLevel FromString(const std::string& level)
{
    if (level == "TRACE") return TRACE;
    if (level == "DEBUG") return DEBUG;
    if (level == "INFO") return INFO;
    if (level == "NOTICE") return NOTICE;
    if (level == "WARN") return WARN;
    if (level == "ERROR") return ERROR;
    if (level == "CRITICAL") return CRITICAL;
    if (level == "FATAL") return FATAL;
    throw std::invalid_argument("invalid log level");
}

class Logger
{
public:
    Logger(std::ostream& os, LogLevel level = INFO)
        : level_(level)
        , os_(os)  // TODO(lhepler) return this to initializer list syntax when gcc-5
        , writer_(&Logger::MessageWriter, this)
    {
#ifdef NDEBUG
        if (level == TRACE)
            throw std::invalid_argument("one cannot simply log TRACE messages in release builds!");
#endif
    }

    Logger(std::ostream& os, const std::string& level) : Logger(os, FromString(level)) {}
    inline LogLevel GetLevel() const { return level_; }
    inline Logger& operator<<(std::unique_ptr<std::ostringstream>&& ptr)
    {
        if (!writer_.joinable()) throw std::runtime_error("this logger is dead!");

        {
            std::lock_guard<std::mutex> g(m_);
            queue_.emplace(std::forward<std::unique_ptr<std::ostringstream>>(ptr));
        }
        pushed_.notify_all();
        return *this;
    }

    inline void Die()
    {
        if (!writer_.joinable()) throw std::runtime_error("this logger is already dead!");

        // place a terminal sentinel for MessageWriter to know it's done
        {
            std::lock_guard<std::mutex> g(m_);
            queue_.emplace(std::unique_ptr<std::ostringstream>());
        }
        pushed_.notify_all();

        // wait for everything to be flushed
        std::unique_lock<std::mutex> lk(m_);
        popped_.wait(lk, [this]() { return queue_.empty(); });
        // endl implicitly flushes, so no need to call os_.flush() here
        // join writer thread
        writer_.join();
    }

    Logger(const Logger& logger) = delete;

    ~Logger()
    {
        if (writer_.joinable()) Die();
    }

    static inline Logger& Default(Logger* logger = nullptr)
    {
        static std::unique_ptr<Logger> logger_(new Logger(std::cerr));
        if (logger) logger_.reset(logger);
        return *logger_;
    }

private:
    void MessageWriter()
    {
        while (true) {
            std::unique_ptr<std::ostringstream> ptr;

            // wait on messages to arrive in the queue_, and pop them off
            {
                std::unique_lock<std::mutex> lk(m_);
                pushed_.wait(lk, [&ptr, this]() {
                    if (queue_.empty()) return false;

                    ptr = std::move(queue_.front());
                    queue_.pop();

                    return true;
                });
            }

            // if we've reached the null terminator, notify flush and stop
            if (!ptr) {
                popped_.notify_all();
                break;
            }

            // otherwise, push the message onto os_
            os_ << ptr->str() << std::endl;

            // and notify flush we delivered a message to os_,
            popped_.notify_all();
        }
    }

private:
    LogLevel level_;
    std::mutex m_;
    std::ostream& os_;
    std::condition_variable popped_;
    std::condition_variable pushed_;
    std::queue<std::unique_ptr<std::ostringstream>> queue_;
    std::thread writer_;
};

class LogMessage
{
public:
    LogMessage(const char* file __attribute__((unused)),
               const char* function __attribute__((unused)),
               unsigned int line __attribute__((unused)), LogLevel level, Logger& logger)
        : logger_(logger)
    {
        using std::chrono::duration_cast;
        using std::chrono::milliseconds;
        using std::chrono::seconds;
        using std::chrono::system_clock;

        if (logger_.GetLevel() > level) return;

        ptr_.reset(new std::ostringstream());

        static const char* delim = " -|- ";

        // get the time, separated into seconds and milliseconds
        const auto now = system_clock::now();
        const auto secs = duration_cast<seconds>(now.time_since_epoch());
        const auto time = system_clock::to_time_t(system_clock::time_point(secs));
        const auto msec = duration_cast<milliseconds>(now.time_since_epoch() - secs).count();

        // format the time and print out the log header to the ostringstream
        // TODO(lhepler) make this std::put_time when we move to gcc-5
        char buf[20];
        std::strftime(buf, 20, "%Y%m%d %T.", std::gmtime(&time));

        (*ptr_) << ">|> " << buf << std::setfill('0') << std::setw(3) << std::to_string(msec)
                << delim << LogLevelRepr(level) << delim << function
#ifndef NDEBUG
                << " at " << file << ':' << line
#endif
                << delim << std::hex << std::showbase << std::this_thread::get_id()
                << std::noshowbase << std::dec << "||" << delim;
    }

    LogMessage(const LogMessage& msg) = delete;

    ~LogMessage()
    {
        if (ptr_) logger_ << std::move(ptr_);
    }

    template <typename T>
    inline LogMessage& operator<<(const T& t)
    {
        if (ptr_) (*ptr_) << t;
        return *this;
    }

private:
    static const char* LogLevelRepr(LogLevel level)
    {
        // by specification these are all of length 10
        switch (level) {
            case TRACE:
                return "TRACE     ";
            case DEBUG:
                return "DEBUG     ";
            case INFO:
                return "INFO      ";
            case NOTICE:
                return "NOTICE    ";
            case WARN:
                return "WARN      ";
            case ERROR:
                return "ERROR     ";
            case CRITICAL:
                return "CRITICAL  ";
            case FATAL:
                return "FATAL     ";
            default:
                return "OTHER     ";
        }
    }

private:
    std::unique_ptr<std::ostringstream> ptr_;
    Logger& logger_;
};

// trace is disabled under Release builds (-DNDEBUG)
#ifdef NDEBUG
#define PBLOGGER_LEVEL(lg, lvl)                                   \
    if (PacBio::CCS::Logging::lvl != PacBio::CCS::Logging::TRACE) \
    PacBio::CCS::Logging::LogMessage(__FILE__, __func__, __LINE__, PacBio::CCS::Logging::lvl, (lg))
#else
#define PBLOGGER_LEVEL(lg, lvl) \
    PacBio::CCS::Logging::LogMessage(__FILE__, __func__, __LINE__, PacBio::CCS::Logging::lvl, (lg))
#endif

#define PBLOGGER_TRACE(lg) PBLOGGER_LEVEL(lg, TRACE)
#define PBLOGGER_DEBUG(lg) PBLOGGER_LEVEL(lg, DEBUG)
#define PBLOGGER_INFO(lg) PBLOGGER_LEVEL(lg, INFO)
#define PBLOGGER_NOTICE(lg) PBLOGGER_LEVEL(lg, NOTICE)
#define PBLOGGER_WARN(lg) PBLOGGER_LEVEL(lg, WARN)
#define PBLOGGER_ERROR(lg) PBLOGGER_LEVEL(lg, ERROR)
#define PBLOGGER_CRITICAL(lg) PBLOGGER_LEVEL(lg, CRITICAL)
#define PBLOGGER_FATAL(lg) PBLOGGER_LEVEL(lg, FATAL)

#define PBLOG_LEVEL(lvl) PBLOGGER_LEVEL(PacBio::CCS::Logging::Logger::Default(), lvl)

#define PBLOG_TRACE PBLOG_LEVEL(TRACE)
#define PBLOG_DEBUG PBLOG_LEVEL(DEBUG)
#define PBLOG_INFO PBLOG_LEVEL(INFO)
#define PBLOG_NOTICE PBLOG_LEVEL(NOTICE)
#define PBLOG_WARN PBLOG_LEVEL(WARN)
#define PBLOG_ERROR PBLOG_LEVEL(ERROR)
#define PBLOG_CRITICAL PBLOG_LEVEL(CRITICAL)
#define PBLOG_FATAL PBLOG_LEVEL(FATAL)

inline void InstallSignalHandlers(Logger& logger = Logger::Default())
{
    using std::raise;
    using std::signal;

    static Logger& logger_ = logger;

    signal(SIGABRT, [](int) {
        {
            PBLOGGER_FATAL(logger_) << "caught SIGABRT";
        }
        logger_.Die();
        signal(SIGABRT, SIG_DFL);
        raise(SIGABRT);
    });
    signal(SIGINT, [](int) {
        {
            PBLOGGER_FATAL(logger_) << "caught SIGINT";
        }
        logger_.Die();
        signal(SIGINT, SIG_DFL);
        raise(SIGINT);
    });
    signal(SIGSEGV, [](int) {
        {
            PBLOGGER_FATAL(logger_) << "caught SIGSEGV";
        }
        logger_.Die();
        signal(SIGSEGV, SIG_DFL);
        raise(SIGSEGV);
    });
    signal(SIGTERM, [](int) {
        {
            PBLOGGER_FATAL(logger_) << "caught SIGTERM";
        }
        logger_.Die();
        signal(SIGTERM, SIG_DFL);
        raise(SIGTERM);
    });
}

}  // namespace Logging
}  // namespace CCS
}  // namespace PacBio
