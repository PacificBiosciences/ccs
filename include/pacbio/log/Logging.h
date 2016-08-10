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
#include <cstdlib>
#include <ctime>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

namespace PacBio {
namespace Logging {

// borrowed a little from primary's SmartEnum
class LogLevel
{
public:
    LogLevel(const unsigned char value) : value_{value} {}
    LogLevel(const std::string& value) : value_{FromString(value)} {}
    enum : unsigned char
    {
        TRACE = 0,
        DEBUG = 1,
        INFO = 2,
        NOTICE = 3,
        WARN = 4,
        ERROR = 5,
        CRITICAL = 6,
        FATAL = 7
    };

    operator unsigned char() const { return value_; }
private:
    inline static LogLevel FromString(const std::string& level)
    {
        if (level == "TRACE") return LogLevel::TRACE;
        if (level == "DEBUG") return LogLevel::DEBUG;
        if (level == "INFO") return LogLevel::INFO;
        if (level == "NOTICE") return LogLevel::NOTICE;
        if (level == "WARN") return LogLevel::WARN;
        if (level == "ERROR") return LogLevel::ERROR;
        if (level == "CRITICAL") return LogLevel::CRITICAL;
        if (level == "FATAL") return LogLevel::FATAL;
        throw std::invalid_argument("invalid log level");
    }

private:
    unsigned char value_;
};

class LoggerConfig : public std::map<LogLevel, std::vector<std::reference_wrapper<std::ostream>>>
{
public:
    LoggerConfig(const std::map<LogLevel, std::vector<std::reference_wrapper<std::ostream>>>& cfg)
        : std::map<LogLevel, std::vector<std::reference_wrapper<std::ostream>>>(cfg)
    {
    }

    LoggerConfig(
        const std::map<std::string, std::vector<std::reference_wrapper<std::ostream>>>& cfg)
    {
        for (const auto& kv : cfg)
            (*this)[kv.first] = kv.second;
    }

    LoggerConfig(std::ostream& os, const LogLevel level = LogLevel::INFO)
    {
        for (size_t i = static_cast<size_t>(level); i < 8; ++i)
            (*this)[static_cast<LogLevel>(i)].push_back(os);
    }

    LoggerConfig(std::ostream& os, const std::string& level) : LoggerConfig(os, LogLevel(level)) {}
};

// necessary fwd decl
class LogMessage;

class Logger
{
public:
    template <typename... Args>
    Logger(Args&&... args)
        : cfg_(std::forward<Args>(args)...), writer_(&Logger::MessageWriter, this)
    {
#ifdef NDEBUG
        if (Handles(LogLevel::TRACE))
            throw std::invalid_argument("one cannot simply log TRACE messages in release builds!");
#endif
    }

    Logger(const Logger& logger) = delete;

    ~Logger()
    {
        if (!writer_.joinable()) throw std::runtime_error("this logger is already dead!");

        // place a terminal sentinel for MessageWriter to know it's done
        {
            std::lock_guard<std::mutex> g(m_);
            queue_.emplace(std::unique_ptr<std::pair<LogLevel, std::ostringstream>>());
        }
        pushed_.notify_all();

        // wait for everything to be flushed
        {
            std::unique_lock<std::mutex> lk(m_);
            popped_.wait(lk, [this]() { return queue_.empty(); });
            // endl implicitly flushes, so no need to call os_.flush() here
        }

        // join writer thread
        writer_.join();
    }

    static inline Logger& Default(Logger* logger = nullptr)
    {
        static std::unique_ptr<Logger> logger_(new Logger(std::cerr));
        if (logger) logger_.reset(logger);
        return *logger_;
    }

private:
    inline bool Handles(const LogLevel level) const
    {
        return cfg_.find(level) != cfg_.end() && !cfg_.at(level).empty();
    }

    inline Logger& operator<<(std::unique_ptr<std::pair<LogLevel, std::ostringstream>>&& ptr)
    {
        if (!writer_.joinable()) throw std::runtime_error("this logger is dead!");

        {
            std::lock_guard<std::mutex> g(m_);
            queue_.emplace(
                std::forward<std::unique_ptr<std::pair<LogLevel, std::ostringstream>>>(ptr));
        }
        pushed_.notify_all();
        return *this;
    }

    void MessageWriter()
    {
        while (true) {
            std::unique_ptr<std::pair<LogLevel, std::ostringstream>> ptr;

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
            const LogLevel level = std::get<0>(*ptr);
            if (cfg_.find(level) != cfg_.end())
                for (const auto& os : cfg_.at(level))
                    os.get() << std::get<1>(*ptr).str() << std::endl;

            // and notify flush we delivered a message to os_,
            popped_.notify_all();
        }
    }

private:
    std::mutex m_;
    LoggerConfig cfg_;
    std::condition_variable popped_;
    std::condition_variable pushed_;
    std::queue<std::unique_ptr<std::pair<LogLevel, std::ostringstream>>> queue_;
    std::thread writer_;

    friend class LogMessage;
};

class LogMessage
{
public:
    LogMessage(const char* file __attribute__((unused)),
               const char* function __attribute__((unused)),
               unsigned int line __attribute__((unused)), const LogLevel level, Logger& logger)
        : logger_(logger)
    {
        using std::chrono::duration_cast;
        using std::chrono::milliseconds;
        using std::chrono::seconds;
        using std::chrono::system_clock;

        if (!logger_.Handles(level)) return;

        ptr_.reset(new std::pair<LogLevel, std::ostringstream>(
            std::piecewise_construct, std::forward_as_tuple(level), std::forward_as_tuple()));

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

        std::get<1>(*ptr_) << ">|> " << buf << std::setfill('0') << std::setw(3)
                           << std::to_string(msec) << delim << LogLevelRepr(level) << delim
                           << function
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
        if (ptr_) std::get<1>(*ptr_) << t;
        return *this;
    }

private:
    static const char* LogLevelRepr(LogLevel level)
    {
        // by specification these are all of length 10
        switch (level) {
            case LogLevel::TRACE:
                return "TRACE     ";
            case LogLevel::DEBUG:
                return "DEBUG     ";
            case LogLevel::INFO:
                return "INFO      ";
            case LogLevel::NOTICE:
                return "NOTICE    ";
            case LogLevel::WARN:
                return "WARN      ";
            case LogLevel::ERROR:
                return "ERROR     ";
            case LogLevel::CRITICAL:
                return "CRITICAL  ";
            case LogLevel::FATAL:
                return "FATAL     ";
            default:
                return "OTHER     ";
        }
    }

private:
    std::unique_ptr<std::pair<LogLevel, std::ostringstream>> ptr_;
    Logger& logger_;
};

// trace is disabled under Release builds (-DNDEBUG)
#ifdef NDEBUG
#define PBLOGGER_LEVEL(lg, lvl)                                   \
    if (PacBio::Logging::lvl != PacBio::Logging::LogLevel::TRACE) \
    PacBio::Logging::LogMessage(__FILE__, __func__, __LINE__, PacBio::Logging::lvl, (lg))
#else
#define PBLOGGER_LEVEL(lg, lvl) \
    PacBio::Logging::LogMessage(__FILE__, __func__, __LINE__, PacBio::Logging::lvl, (lg))
#endif

#define PBLOGGER_TRACE(lg) PBLOGGER_LEVEL(lg, LogLevel::TRACE)
#define PBLOGGER_DEBUG(lg) PBLOGGER_LEVEL(lg, LogLevel::DEBUG)
#define PBLOGGER_INFO(lg) PBLOGGER_LEVEL(lg, LogLevel::INFO)
#define PBLOGGER_NOTICE(lg) PBLOGGER_LEVEL(lg, LogLevel::NOTICE)
#define PBLOGGER_WARN(lg) PBLOGGER_LEVEL(lg, LogLevel::WARN)
#define PBLOGGER_ERROR(lg) PBLOGGER_LEVEL(lg, LogLevel::ERROR)
#define PBLOGGER_CRITICAL(lg) PBLOGGER_LEVEL(lg, LogLevel::CRITICAL)
#define PBLOGGER_FATAL(lg) PBLOGGER_LEVEL(lg, LogLevel::FATAL)

#define PBLOG_LEVEL(lvl) PBLOGGER_LEVEL(PacBio::Logging::Logger::Default(), lvl)

#define PBLOG_TRACE PBLOG_LEVEL(LogLevel::TRACE)
#define PBLOG_DEBUG PBLOG_LEVEL(LogLevel::DEBUG)
#define PBLOG_INFO PBLOG_LEVEL(LogLevel::INFO)
#define PBLOG_NOTICE PBLOG_LEVEL(LogLevel::NOTICE)
#define PBLOG_WARN PBLOG_LEVEL(LogLevel::WARN)
#define PBLOG_ERROR PBLOG_LEVEL(LogLevel::ERROR)
#define PBLOG_CRITICAL PBLOG_LEVEL(LogLevel::CRITICAL)
#define PBLOG_FATAL PBLOG_LEVEL(LogLevel::FATAL)

inline void InstallSignalHandlers(Logger& logger = Logger::Default())
{
    using std::raise;
    using std::signal;

    static Logger& logger_ = logger;

    std::set_terminate([]() {
        if (auto eptr = std::current_exception()) {
            try {
                std::rethrow_exception(eptr);
            } catch (const std::exception& e) {
                PBLOGGER_FATAL(logger_) << "caught exception: \"" << e.what() << '"';
            } catch (...) {
                PBLOGGER_FATAL(logger_) << "caught unknown exception type";
            }
        }
        // call the SIGABRT handler (below)
        std::abort();
    });
    signal(SIGABRT, [](int) {
        {
            PBLOGGER_FATAL(logger_) << "caught SIGABRT";
        }
        logger_.~Logger();
        signal(SIGABRT, SIG_DFL);
        raise(SIGABRT);
    });
    signal(SIGINT, [](int) {
        {
            PBLOGGER_FATAL(logger_) << "caught SIGINT";
        }
        logger_.~Logger();
        signal(SIGINT, SIG_DFL);
        raise(SIGINT);
    });
    signal(SIGSEGV, [](int) {
        {
            PBLOGGER_FATAL(logger_) << "caught SIGSEGV";
        }
        logger_.~Logger();
        signal(SIGSEGV, SIG_DFL);
        raise(SIGSEGV);
    });
    signal(SIGTERM, [](int) {
        {
            PBLOGGER_FATAL(logger_) << "caught SIGTERM";
        }
        logger_.~Logger();
        signal(SIGTERM, SIG_DFL);
        raise(SIGTERM);
    });
}

}  // namespace Logging
}  // namespace PacBio
