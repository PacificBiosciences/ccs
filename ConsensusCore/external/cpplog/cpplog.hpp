#pragma once

#ifndef _CPPLOG_H
#define _CPPLOG_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <ctime>
#include <vector>
#include <cstdlib>
#include <streambuf>
#include <ostream>

// The following #define's will change the behaviour of this library.
//      #define CPPLOG_FILTER_LEVEL     <level>
//          Prevents all log messages with level less than <level> from being emitted.
//
//      #define CPPLOG_SYSTEM_IDS
//          Enables capturing of the Process and Thread ID.
//
//      #define CPPLOG_USE_SYSCALL_FOR_THREAD_ID
//          Tries to use syscall(SYS_gettid) to get the current Thread ID.  This is
//          especially useful on systems where boost's get_current_thread_id()
//          will return a structure (usually pthread_t), which isn't incredibly
//          useful.
//          NOTE: Only useful if you also #define CPPLOG_SYSTEM_IDS
//
//      #define CPPLOG_THREADING
//          Enables threading (BackgroundLogger).  Note that defining this or
//          CPPLOG_SYSTEM_IDS introduces a dependency on Boost;
//          this means that the library is no longer truly header-only.
//
//      #define CPPLOG_HELPER_MACROS
//          Enables inclusion of the CHECK_* macros.
//
//      #define CPPLOG_FATAL_EXIT
//          Causes a fatal error to exit() the process.
//
//      #define CPPLOG_FATAL_EXIT_DEBUG
//          Causes a fatal error to exit() the process if in debug mode.
//
//      # define CPPLOG_USE_OLD_BOOST
//          Use the old Boost namespace for interprocess::ipcdetail.  Define
//          this if you're using version 1.47 of Boost or earlier.

// ------------------------------- DEFINITIONS -------------------------------

// Severity levels:
// Note: Giving a value for CPPLOG_FILTER_LEVEL will log all messages at
//       or above that level.
//  0 - Trace
//  1 - Debug
//  2 - Info
//  3 - Warning
//  4 - Error
//  5 - Fatal (always logged)

#define LL_TRACE    0
#define LL_DEBUG    1
#define LL_INFO     2
#define LL_WARN     3
#define LL_ERROR    4
#define LL_FATAL    5


// ------------------------------ CONFIGURATION ------------------------------

//#define CPPLOG_FILTER_LEVEL               LL_WARN
//#define CPPLOG_SYSTEM_IDS
//#define CPPLOG_USE_SYSCALL_FOR_THREAD_ID
//#define CPPLOG_THREADING
#define CPPLOG_HELPER_MACROS
#define CPPLOG_FATAL_EXIT
//#define CPPLOG_FATAL_EXIT_DEBUG
//#define CPPLOG_USE_OLD_BOOST


// ---------------------------------- CODE -----------------------------------

#ifdef CPPLOG_SYSTEM_IDS
#include <boost/interprocess/detail/os_thread_functions.hpp>
#ifdef CPPLOG_USE_SYSCALL_FOR_THREAD_ID
#include <unistd.h>
#include <sys/syscall.h>
#endif
#endif

#ifdef CPPLOG_THREADING
#include <boost/thread.hpp>
#include "concurrent_queue.hpp"
#endif

#ifdef _WIN32
#include "outputdebugstream.hpp"
#endif

#ifdef CPPLOG_WITH_SCRIBE_LOGGER
#include "scribestream.hpp"
#endif

// If we don't have a level defined, set it to CPPLOG_LEVEL_DEBUG (log all except trace statements)
#ifndef CPPLOG_FILTER_LEVEL
#define CPPLOG_FILTER_LEVEL LL_DEBUG
#endif


// The general concept for how logging works:
//  - Every call to LOG(LEVEL, logger) works as follows:
//      - Instantiates an object of type LogMessage.
//      - LogMessage's constructor captures __FILE__, __LINE__, severity and our output logger.
//      - LogMessage exposes a function getStream(), which is an ostringstream-style stream that
//        client code can write debug information into.
//      - When the sendToLogger() method of a LogMessage is called, all the buffered data in the
//        messages' stream is sent to the specified logger.
//      - When a LogMessage's destructor is called, it calls sendToLogger() to send all
//        remaining data.


namespace cpplog
{
    // Our log level type.
    // NOTE: When C++11 becomes widely supported, convert this to "enum class LogLevel".
    typedef unsigned int loglevel_t;

    // Helper functions.  Stuck these in their own namespace for simplicity.
    namespace helpers
    {
        // Gets the filename from a path.
        inline static const char* fileNameFromPath(const char* filePath)
        {
            const char* fileName = strrchr(filePath, '/');
#if defined(_WIN32)
            if( !fileName )
                fileName = strrchr(filePath, '\\');
#endif
            return fileName ? fileName + 1 : filePath;
        }

        // Thread-safe version of localtime()
        inline bool slocaltime(::tm* const out, const ::time_t* const in)
        {
#if defined(_WIN32) && defined(_MSC_VER)
            return ::localtime_s(out, in) == 0;
#elif defined(__MINGW32__)
            // Warning - not entirely thread safe on MinGW
            ::tm * localOut = ::localtime(in);
            if( localOut )
            {
                ::memcpy(out, localOut, sizeof(::tm));
            }
            return localOut != NULL;
#else
            // Default to SUSv2 (libc >= 5.2.5) function.
            return ::localtime_r(in, out) != NULL;
#endif
        }

        // Thread-safe version of gmtime()
        inline bool sgmtime(::tm* const out, const ::time_t* const in)
        {
#if defined(_WIN32) && defined(_MSC_VER)
            return ::gmtime_s(out, in) == 0;
#elif defined(__MINGW32__)
            // Warning - not entirely thread safe on MinGW
            ::tm * localOut = ::gmtime(in);
            if( localOut )
            {
                ::memcpy(out, localOut, sizeof(::tm));
            }
            return localOut != NULL;
#else
            // Default to SUSv2 (libc >= 5.2.5) function.
            return ::gmtime_r(in, out) != NULL;
#endif
        }

        // Below we have a bunch of macros, typedefs and such that make getting our
        // current process/thread ID simpler.
#ifdef CPPLOG_SYSTEM_IDS
#ifdef CPPLOG_USE_OLD_BOOST
        typedef boost::interprocess::detail::OS_process_id_t    process_id_t;

        inline process_id_t get_process_id()
        {
            return boost::interprocess::detail::get_current_process_id();
        }
#else
        typedef boost::interprocess::ipcdetail::OS_process_id_t process_id_t;

        inline process_id_t get_process_id()
        {
            return boost::interprocess::ipcdetail::get_current_process_id();
        }
#endif

#ifdef CPPLOG_USE_SYSCALL_FOR_THREAD_ID
        typedef unsigned long                                   thread_id_t;

        inline thread_id_t get_thread_id()
        {
            return static_cast<unsigned long>(syscall(SYS_gettid));
        }

        inline void print_thread_id(std::ostream& stream, thread_id_t thread_id)
        {
            stream << std::setfill('0') << std::setw(8) << std::hex
                   << thread_id;
        }
#else   // CPPLOG_USE_SYSCALL_FOR_THREAD_ID
#ifdef CPPLOG_USE_OLD_BOOST
        typedef boost::interprocess::detail::OS_thread_id_t     thread_id_t;

        inline thread_id_t get_thread_id()
        {
            return boost::interprocess::detail::get_current_thread_id();
        }
#else
        typedef boost::interprocess::ipcdetail::OS_thread_id_t  thread_id_t;

        inline thread_id_t get_thread_id()
        {
            return boost::interprocess::ipcdetail::get_current_thread_id();
        }
#endif  // CPPLOG_USE_OLD_BOOST

        // This function lets us print a thread ID in all cases, including on
        // platforms where it's actually a structure (pthread_t, I'm looking
        // at you...).  Note that this kinda-sorta assumes a little-endian
        // architecture, if we want meaningful results.  Not super important,
        // though, since the address of a structure isn't actually that useful.
        // TODO: I might rewrite this using templates to print properly if it's
        // an unsigned long, and fall back to this implementation otherwise.
        inline void print_thread_id(std::ostream& stream, thread_id_t thread_id)
        {
            unsigned char* sptr = static_cast<unsigned char*>(
                                    static_cast<void*>(&thread_id)
                                  );
            for( size_t i = sizeof(thread_id_t); i != 0; i-- )
            {
                stream << std::setfill('0') << std::setw(2) << std::hex
                       << static_cast<unsigned>(sptr[i - 1]);
            }
        }
#endif  // CPPLOG_USE_SYSCALL_FOR_THREAD_ID
#endif  // CPPLOG_SYSTEM_IDS

        // Simple class that allows us to evaluate a stream to void - prevents compiler errors.
        class VoidStreamClass
        {
        public:
            VoidStreamClass() { }
            void operator&(std::ostream&) { }
        };

        // fixed_streambuf is a minimal implementation around std::basic_streambuf
        // with a fixed size backing buffer. It implements additional functionality
        // needed by cpplog and exposes the backing buffer in a safe way via c_str().
        // This makes it possible to avoid extra copying.
        class fixed_streambuf : public std::basic_streambuf<char, std::char_traits<char> >
        {
        private:
            // Constant.
            static const size_t k_logBufferCapacity = 20000;
            // Leave room for terminating null character in case buffer fills up.
            char m_buffer[k_logBufferCapacity+1];

        public:
            fixed_streambuf()
            {
                // Use allocated buffer as backing store.
                setp(m_buffer, m_buffer + k_logBufferCapacity);
                // Insert terminator at buffer end.
                m_buffer[k_logBufferCapacity] = '\0';
            }

            std::streamsize length()   const { return pptr() - pbase();       }
            std::streamsize capacity() const { return k_logBufferCapacity;    }
            bool empty()               const { return length() == 0;          }
            bool full()                const { return length() == capacity(); }

            // Unput one character.
            int_type sunputc()
            {
                if( (!pptr()) || (pptr() == pbase()) )
                    return pbackfail();

                pbump(-1);

                // This is safe because *epptr() always is '\0' and inside
                // the backing buffer.
                return traits_type::to_int_type(*(pptr()+1));
            }

            // Peek at last inserted character.
            int peek() const
            {
                if( (!pptr()) || (pptr() == pbase()) )
                    return std::char_traits<char>::eof();

                return static_cast<int>(*(pptr()-1));
            }

            const char* c_str() const
            {
                // Add terminating null character.
                // This is safe even if the buffer is full to its capacity since
                // epptr() is inside the backing buffer.
                *pptr() = '\0';
                return pbase();
            }
        };
    }

    // Logger data.  This is sent to a logger when a LogMessage is Flush()'ed, or
    // when the destructor is called.
    struct LogData
    {

        // Our streambuf & stream to log data to.
        helpers::fixed_streambuf streamBuffer;
        std::ostream stream;

        // Captured data.
        unsigned int level;
        unsigned long line;
        const char* fullPath;
        const char* fileName;
        time_t messageTime;
        ::tm utcTime;

#ifdef CPPLOG_SYSTEM_IDS
        // Process/thread ID.
        helpers::process_id_t processId;
        helpers::thread_id_t  threadId;
#endif


        // Constructor that initializes our stream.
        LogData(loglevel_t logLevel)
            : streamBuffer(), stream(&streamBuffer), level(logLevel)
#ifdef CPPLOG_SYSTEM_IDS
              , processId(0), threadId(0)
#endif
        {
        }

        virtual ~LogData()
        { }
    };

    // Base interface for a logger.
    class BaseLogger
    {
    public:
        // All loggers must provide an interface to log a message to.
        // The return value of this function indicates whether to delete
        // the log message.
        virtual bool sendLogMessage(LogData* logData) = 0;

        virtual ~BaseLogger() { }
    };

    // Log message - this is instantiated upon every call to LOG(logger)
    class LogMessage
    {
    private:
        BaseLogger*     m_logger;
        bool            m_flushed;
        bool            m_deleteMessage;

    protected:
        LogData*        m_logData;

    private:
        // Flag for if a fatal message has been logged already.
        // This prevents us from calling exit(), which calls something,
        // which then logs a fatal message, which cause an infinite loop.
        // TODO: this should probably be thread-safe...
        static bool getSetFatal(bool get, bool val)
        {
            static bool m_fatalFlag = false;

            if( !get )
                m_fatalFlag = val;

            return m_fatalFlag;
        }

    public:
        LogMessage(const char* file, unsigned int line, loglevel_t logLevel, BaseLogger* outputLogger, bool useDefaultLogFormat=true)
            : m_logger(outputLogger)
        {
            Init(file, line, logLevel, useDefaultLogFormat);
        }

        LogMessage(const char* file, unsigned int line, loglevel_t logLevel, BaseLogger& outputLogger, bool useDefaultLogFormat=true)
            : m_logger(&outputLogger)
        {
            Init(file, line, logLevel, useDefaultLogFormat);
        }

        virtual ~LogMessage()
        {
            Flush();

            if( m_deleteMessage )
            {
                delete m_logData;
            }
        }

        inline std::ostream& getStream()
        {
            return m_logData->stream;
        }

    protected:
        virtual void InitLogMessage()
        {
            // Log process ID and thread ID.
#ifdef CPPLOG_SYSTEM_IDS
            m_logData->stream << "["
                        << std::right << std::setfill('0') << std::setw(8) << std::hex
                        << m_logData->processId << ".";
            helpers::print_thread_id(m_logData->stream, m_logData->threadId);
            m_logData->stream << "] ";
#endif

            m_logData->stream << std::setfill(' ') << std::setw(5) << std::left << std::dec
                        << LogMessage::getLevelName(m_logData->level) << " - "
                        << m_logData->fileName << "(" << m_logData->line << "): ";
        }

    private:
        void Init(const char* file, unsigned int line, loglevel_t logLevel, bool useDefaultLogFormat=true)
        {
            m_logData = new LogData(logLevel);
            m_flushed = false;
            m_deleteMessage = false;

            // Capture data.
            m_logData->fullPath     = file;
            m_logData->fileName     = cpplog::helpers::fileNameFromPath(file);
            m_logData->line         = line;
            m_logData->messageTime  = ::time(NULL);

            // Get current time.
            ::tm gmt;
            cpplog::helpers::sgmtime(&gmt, &m_logData->messageTime);
            memcpy(&m_logData->utcTime, &gmt, sizeof(tm));

#ifdef CPPLOG_SYSTEM_IDS
            // Get process/thread ID.
            m_logData->processId    = helpers::get_process_id();
            m_logData->threadId     = helpers::get_thread_id();
#endif // CPPLOG_SYSTEM_IDS

            if( useDefaultLogFormat )
            {
                InitLogMessage();
            }
        }

        void Flush()
        {
            if( !m_flushed )
            {
                // Insert newline, if needed.
                helpers::fixed_streambuf* const sb = &m_logData->streamBuffer;
                if( sb->peek() != '\n' )
                {
                    // If buffer is full, remove last char to leave room for newline.
                    if( sb->full() )
                        sb->sunputc();

                    // Insert newline
                    sb->sputc('\n');
                }

                // Save the log level.
                loglevel_t savedLogLevel = m_logData->level;

                // Send the message, set flushed=true.
                m_deleteMessage = m_logger->sendLogMessage(m_logData);
                m_flushed = true;

                // Note: We cannot touch m_logData after the above call.  By the
                // time it returns, we have to assume it has already been freed.

                // If this is a fatal message...
                if( savedLogLevel == LL_FATAL && !getSetFatal(true, true) )
                {
                    // Set our fatal flag.
                    getSetFatal(false, true);

#ifdef _DEBUG
// Only exit in debug mode if CPPLOG_FATAL_EXIT_DEBUG is set.
#if defined(CPPLOG_FATAL_EXIT_DEBUG) || defined(CPPLOG_FATAL_EXIT)
                    std::exit(1);
#endif
#else //!_DEBUG
#ifdef CPPLOG_FATAL_EXIT_DEBUG
                    std::exit(1)
#endif
#endif
                }
            }
        }

    public:
        static const char* getLevelName(loglevel_t level)
        {
            switch( level )
            {
                case LL_TRACE:
                    return "TRACE";
                case LL_DEBUG:
                    return "DEBUG";
                case LL_INFO:
                    return "INFO";
                case LL_WARN:
                    return "WARN";
                case LL_ERROR:
                    return "ERROR";
                case LL_FATAL:
                    return "FATAL";
                default:
                    return "OTHER";
            };
        };
    };

    // Generic class - logs to a given std::ostream.
    class OstreamLogger : public BaseLogger
    {
    protected:
        std::ostream&   m_logStream;

    public:
        OstreamLogger(std::ostream& outStream)
            : m_logStream(outStream)
        { }

        virtual bool sendLogMessage(LogData* logData)
        {
            helpers::fixed_streambuf* const sb = &logData->streamBuffer;
            m_logStream.write(sb->c_str(), sb->length());
            m_logStream << std::flush;

            return true;
        }

        virtual ~OstreamLogger() { }
    };

    // Simple implementation - logs to stderr.
    class StdErrLogger : public OstreamLogger
    {
    public:
        StdErrLogger()
            : OstreamLogger(std::cerr)
        { }
    };

    // Simple implementation - logs to a string, provides the ability to get that string.
    class StringLogger : public OstreamLogger
    {
    private:
        std::ostringstream  m_stream;
    public:
        StringLogger()
            : OstreamLogger(m_stream)
        { }

        std::string getString()
        {
            return m_stream.str();
        }

        void clear()
        {
            m_stream.str("");
            m_stream.clear();
        }
    };

#ifdef _WIN32
    class OutputDebugStringLogger : public OstreamLogger
    {
    private:
        dbgwin_stream m_stream;
    public:
        OutputDebugStringLogger() : OstreamLogger(m_stream)
        { }
    };
#endif

    // Log to file.
    class FileLogger : public OstreamLogger
    {
    private:
        std::string     m_path;
        std::ofstream   m_outStream;

    public:
        FileLogger(std::string logFilePath)
            : OstreamLogger(m_outStream), m_path(logFilePath), m_outStream(logFilePath.c_str(), std::ios_base::out)
        {
        }

        FileLogger(std::string logFilePath, bool append)
            : OstreamLogger(m_outStream), m_path(logFilePath), m_outStream(logFilePath.c_str(), append ? std::ios_base::app : std::ios_base::out)
        {
        }
    };

    // Log to file, rotate when the log reaches a given size.
    class SizeRotateFileLogger : public OstreamLogger
    {
    public:
        typedef void (*pfBuildFileName)(unsigned long logNumber, std::string& newFileName, void* context);

    private:
        std::streamoff  m_maxSize;
        unsigned long   m_logNumber;

        SizeRotateFileLogger::pfBuildFileName m_buildFunc;
        void*           m_context;

        std::ofstream   m_outStream;

    public:
        SizeRotateFileLogger(pfBuildFileName nameFunc, std::streamoff maxSize)
            : OstreamLogger(m_outStream), m_maxSize(maxSize), m_logNumber(0),
              m_buildFunc(nameFunc), m_context(NULL),
              m_outStream()
        {
            // "Rotate" to open our initial log.
            RotateLog();
        }

        SizeRotateFileLogger(pfBuildFileName nameFunc, void* context,
                std::streamoff maxSize)
            : OstreamLogger(m_outStream), m_maxSize(maxSize), m_logNumber(0),
              m_buildFunc(nameFunc), m_context(context),
              m_outStream()
        {
            // "Rotate" to open our initial log.
            RotateLog();
        }

        virtual ~SizeRotateFileLogger()
        { }

        virtual bool sendLogMessage(LogData* logData)
        {
            // Call the actual logger.
            bool deleteMessage = OstreamLogger::sendLogMessage(logData);

            // Check if we're over our limit.
            if( m_outStream.tellp() > m_maxSize )
            {
                // Yep, increment our log number and rotate.
                m_logNumber++;
                m_outStream << std::flush;

                RotateLog();
            }

            return deleteMessage;
        }


    private:
        void RotateLog()
        {
            // Build the file name.
            std::string newFileName;
            m_buildFunc(m_logNumber, newFileName, m_context);

            // Close old file, open new file.
            m_outStream.close();
            m_outStream.open(newFileName.c_str(), std::ios_base::out);
        }
    };

    // Log to file, rotate every "x" seconds.
    class TimeRotateFileLogger : public OstreamLogger
    {
    public:
        typedef void (*pfBuildFileName)(::tm* time, unsigned long logNumber,
                                        std::string& newFileName, void* context);

    private:
        unsigned long   m_rotateInterval;
        ::time_t        m_lastRotateTime;
        unsigned long   m_logNumber;

        cpplog::TimeRotateFileLogger::pfBuildFileName m_buildFunc;
        void* m_context;

        std::ofstream   m_outStream;

    public:
        TimeRotateFileLogger(pfBuildFileName nameFunc, unsigned long intervalSeconds)
            : OstreamLogger(m_outStream), m_rotateInterval(intervalSeconds), m_logNumber(0),
              m_buildFunc(nameFunc), m_context(NULL)
        {
            // "Rotate" to open our initial log.
            RotateLog(::time(NULL));
        }

        TimeRotateFileLogger(pfBuildFileName nameFunc, void* context, unsigned long intervalSeconds)
            : OstreamLogger(m_outStream), m_rotateInterval(intervalSeconds), m_logNumber(0),
              m_buildFunc(nameFunc), m_context(context)
        {
            // "Rotate" to open our initial log.
            RotateLog(::time(NULL));
        }

        virtual ~TimeRotateFileLogger()
        {
        }

        virtual bool sendLogMessage(LogData* logData)
        {
            // Get the current time.
            ::time_t currTime;
            ::time(&currTime);

            unsigned long timeDiff = static_cast<unsigned long>(
                                        difftime(currTime, m_lastRotateTime)
                                     );

            // Is the difference greater than our number of seconds?
            if( timeDiff > m_rotateInterval )
            {
                // Yep, increment our log number and rotate.
                m_logNumber++;
                m_outStream << std::flush;

                RotateLog(currTime);
            }

            // Call the actual logger.
            return OstreamLogger::sendLogMessage(logData);
        }

    private:
        void RotateLog(time_t currTime)
        {
            // Get the current time.
            ::tm timeInfo;
            cpplog::helpers::slocaltime(&timeInfo, &currTime);

            // Build a new file name.
            std::string newFileName;
            m_buildFunc(&timeInfo, m_logNumber, newFileName, m_context);

            // Close old file, open new file.
            m_outStream.close();
            m_outStream.open(newFileName.c_str(), std::ios_base::out);

            // Reset the rotate time.
            ::time(&m_lastRotateTime);
        }
    };

#ifdef CPPLOG_WITH_SCRIBE_LOGGER
    // Given a Scribe node, will send log messages there with the given category.
    class ScribeLogger : public OstreamLogger
    {
    private:
        scribe_stream m_outStream;
    public:
        ScribeLogger(std::string host, unsigned short port, std::string category, int timeout)
            : OstreamLogger(m_outStream)
        {
            m_outStream.open(host, port, category, timeout);
        }
    };
#endif

    // Tee logger - given two loggers, will forward a message to both.
    class TeeLogger : public BaseLogger
    {
    private:
        BaseLogger* m_logger1;
        BaseLogger* m_logger2;

        bool        m_logger1Owned;
        bool        m_logger2Owned;

    public:
        TeeLogger(BaseLogger* one, BaseLogger* two)
            : m_logger1(one), m_logger2(two),
              m_logger1Owned(false), m_logger2Owned(false)
        { }

        TeeLogger(BaseLogger* one, bool ownOne, BaseLogger* two, bool ownTwo)
            : m_logger1(one), m_logger2(two),
              m_logger1Owned(ownOne), m_logger2Owned(ownTwo)
        { }

        TeeLogger(BaseLogger& one, BaseLogger& two)
            : m_logger1(&one), m_logger2(&two),
              m_logger1Owned(false), m_logger2Owned(false)
        { }

        TeeLogger(BaseLogger& one, bool ownOne, BaseLogger& two, bool ownTwo)
            : m_logger1(&one), m_logger2(&two),
              m_logger1Owned(ownOne), m_logger2Owned(ownTwo)
        { }

        ~TeeLogger()
        {
            if( m_logger1Owned )
                delete m_logger1;
            if( m_logger2Owned )
                delete m_logger2;
        }

        virtual bool sendLogMessage(LogData* logData)
        {
            bool deleteMessage = true;

            deleteMessage = deleteMessage && m_logger1->sendLogMessage(logData);
            deleteMessage = deleteMessage && m_logger2->sendLogMessage(logData);

            return deleteMessage;
        }
    };

    // Multiplex logger - will forward a log message to all loggers.
    class MultiplexLogger : public BaseLogger
    {
        struct LoggerInfo
        {
            BaseLogger* logger;
            bool        owned;

            LoggerInfo(BaseLogger* l, bool o)
                : logger(l), owned(o)
            { }
        };
        std::vector<LoggerInfo> m_loggers;

    public:
        MultiplexLogger()
        { }

        MultiplexLogger(BaseLogger* one)
        {
            m_loggers.push_back(LoggerInfo(one, false));
        }

        MultiplexLogger(BaseLogger& one)
        {
            m_loggers.push_back(LoggerInfo(&one, false));
        }

        MultiplexLogger(BaseLogger* one, bool owned)
        {
            m_loggers.push_back(LoggerInfo(one, owned));
        }

        MultiplexLogger(BaseLogger& one, bool owned)
        {
            m_loggers.push_back(LoggerInfo(&one, owned));
        }

        MultiplexLogger(BaseLogger* one, BaseLogger* two)
        {
            m_loggers.push_back(LoggerInfo(one, false));
            m_loggers.push_back(LoggerInfo(two, false));
        }

        MultiplexLogger(BaseLogger* one, bool ownOne, BaseLogger* two, bool ownTwo)
        {
            m_loggers.push_back(LoggerInfo(one, ownOne));
            m_loggers.push_back(LoggerInfo(two, ownTwo));
        }

        MultiplexLogger(BaseLogger& one, bool ownOne, BaseLogger& two, bool ownTwo)
        {
            m_loggers.push_back(LoggerInfo(&one, ownOne));
            m_loggers.push_back(LoggerInfo(&two, ownTwo));
        }

        ~MultiplexLogger()
        {
            for( std::vector<LoggerInfo>::iterator It = m_loggers.begin();
                 It != m_loggers.end();
                 It++ )
            {
                if( (*It).owned )
                    delete (*It).logger;
            }
        }

        void addLogger(BaseLogger* logger)      { m_loggers.push_back(LoggerInfo(logger, false)); }
        void addLogger(BaseLogger& logger)      { m_loggers.push_back(LoggerInfo(&logger, false)); }

        void addLogger(BaseLogger* logger, bool owned)      { m_loggers.push_back(LoggerInfo(logger, owned)); }
        void addLogger(BaseLogger& logger, bool owned)      { m_loggers.push_back(LoggerInfo(&logger, owned)); }

        virtual bool sendLogMessage(LogData* logData)
        {
            bool deleteMessage = true;

            for( std::vector<LoggerInfo>::iterator It = m_loggers.begin();
                 It != m_loggers.end();
                 It++ )
            {
                deleteMessage = deleteMessage && (*It).logger->sendLogMessage(logData);
            }

            return deleteMessage;
        }
    };

    // Filtering logger.  Will not forward all messages less than a given level.
    class FilteringLogger : public BaseLogger
    {
    private:
        loglevel_t      m_lowestLevelAllowed;
        BaseLogger*     m_forwardTo;
        bool            m_owned;

    public:
        FilteringLogger(loglevel_t level, BaseLogger* forwardTo)
            : m_lowestLevelAllowed(level), m_forwardTo(forwardTo), m_owned(false)
        { }

        FilteringLogger(loglevel_t level, BaseLogger& forwardTo)
            : m_lowestLevelAllowed(level), m_forwardTo(&forwardTo), m_owned(false)
        { }

        FilteringLogger(loglevel_t level, BaseLogger* forwardTo, bool owned)
            : m_lowestLevelAllowed(level), m_forwardTo(forwardTo), m_owned(owned)
        { }

        FilteringLogger(loglevel_t level, BaseLogger& forwardTo, bool owned)
            : m_lowestLevelAllowed(level), m_forwardTo(&forwardTo), m_owned(owned)
        { }

        ~FilteringLogger()
        {
            if( m_owned )
                delete m_forwardTo;
        }

        void SetLevel(loglevel_t allowed)
        {
            m_lowestLevelAllowed = allowed;
        }

        virtual bool sendLogMessage(LogData* logData)
        {
            if( logData->level >= m_lowestLevelAllowed )
                return m_forwardTo->sendLogMessage(logData);
            else
                return true;
        }
    };

    // Logger that moves all processing of log messages to a background thread.
    // Only include if we have support for threading.
#ifdef CPPLOG_THREADING
    class BackgroundLogger : public BaseLogger
    {
    private:
        BaseLogger*                 m_forwardTo;
        concurrent_queue<LogData*>  m_queue;

        boost::thread               m_backgroundThread;
        LogData*                    m_dummyItem;

        void backgroundFunction()
        {
            LogData* nextLogEntry;
            bool deleteMessage = true;

            do
            {
                m_queue.wait_and_pop(nextLogEntry);

                if( nextLogEntry != m_dummyItem )
                    deleteMessage = m_forwardTo->sendLogMessage(nextLogEntry);

                if( deleteMessage )
                    delete nextLogEntry;
            } while( nextLogEntry != m_dummyItem );
        }

        void Init()
        {
            // Create dummy item.
            m_dummyItem = new LogData(LL_TRACE);

            // And create background thread.
            m_backgroundThread = boost::thread(&BackgroundLogger::backgroundFunction, this);
        }

    public:
        BackgroundLogger(BaseLogger* forwardTo)
            : m_forwardTo(forwardTo)
        {
            Init();
        }

        BackgroundLogger(BaseLogger& forwardTo)
            : m_forwardTo(&forwardTo)
        {
            Init();
        }

        void Stop()
        {
            // Push our "dummy" item on the queue ...
            m_queue.push(m_dummyItem);

            // ... and wait for thread to terminate.
            m_backgroundThread.join();

            // NOTE: The loop will free the dummy item for us, we can ignore it.
        }

        ~BackgroundLogger()
        {
            Stop();
        }

        virtual bool sendLogMessage(LogData* logData)
        {
            m_queue.push(logData);

            // Don't delete - the background thread should handle this.
            return false;
        }

    };

#endif

    // Seperate namespace for loggers that use templates.
    namespace templated
    {
        // Filtering logger that accepts the level as a template parameter.
        // This will be slightly faster at runtime, as the if statement will
        // be performed on a constant value, as opposed to needing a memory
        // lookup (as with FilteringLogger)
        template <loglevel_t lowestLevel = LL_TRACE>
        class TFilteringLogger : public BaseLogger
        {
            BaseLogger* m_forwardTo;

        public:
            TFilteringLogger(BaseLogger* forwardTo)
                : m_forwardTo(forwardTo)
            { }

            virtual bool sendLogMessage(LogData* logData)
            {
                if( logData->level >= lowestLevel )
                    return m_forwardTo->sendLogMessage(logData);
                else
                    return true;
            }
        };

        // TODO: Implement others?
    }
}

// Our logging macros.

// Default macros - log, and don't log something.
// Allow custom log message formatting
#ifndef LOG_LEVEL
#define LOG_LEVEL(level, logger)    cpplog::LogMessage(__FILE__, __LINE__, (level), logger).getStream()
#endif
#define LOG_NOTHING(level, logger)  true ? (void)0 : cpplog::helpers::VoidStreamClass() & LOG_LEVEL(level, logger)

// Series of debug macros, depending on what we log.
#if CPPLOG_FILTER_LEVEL <= LL_TRACE
#define LOG_TRACE(logger)   LOG_LEVEL(LL_TRACE, logger)
#else
#define LOG_TRACE(logger)   LOG_NOTHING(LL_TRACE, logger)
#endif

#if CPPLOG_FILTER_LEVEL <= LL_DEBUG
#define LOG_DEBUG(logger)   LOG_LEVEL(LL_DEBUG, logger)
#else
#define LOG_DEBUG(logger)   LOG_NOTHING(LL_DEBUG, logger)
#endif

#if CPPLOG_FILTER_LEVEL <= LL_INFO
#define LOG_INFO(logger)    LOG_LEVEL(LL_INFO, logger)
#else
#define LOG_INFO(logger)    LOG_NOTHING(LL_INFO, logger)
#endif

#if CPPLOG_FILTER_LEVEL <= LL_WARN
#define LOG_WARN(logger)    LOG_LEVEL(LL_WARN, logger)
#else
#define LOG_WARN(logger)    LOG_NOTHING(LL_WARN, logger)
#endif

#if CPPLOG_FILTER_LEVEL <= LL_ERROR
#define LOG_ERROR(logger)   LOG_LEVEL(LL_ERROR, logger)
#else
#define LOG_ERROR(logger)   LOG_NOTHING(LL_ERROR, logger)
#endif

// Note: Always logged.
#define LOG_FATAL(logger)   LOG_LEVEL(LL_FATAL, logger)



// Debug macros - only logged in debug mode.
#ifdef _DEBUG
#define DLOG_TRACE(logger)  LOG_TRACE(logger)
#define DLOG_DEBUG(logger)  LOG_DEBUG(logger)
#define DLOG_INFO(logger)   LOG_INFO(logger)
#define DLOG_WARN(logger)   LOG_WARN(logger)
#define DLOG_ERROR(logger)  LOG_ERROR(logger)
#else
#define DLOG_TRACE(logger)  LOG_NOTHING(LL_TRACE, logger)
#define DLOG_DEBUG(logger)  LOG_NOTHING(LL_DEBUG, logger)
#define DLOG_INFO(logger)   LOG_NOTHING(LL_INFO,  logger)
#define DLOG_WARN(logger)   LOG_NOTHING(LL_WARN,  logger)
#define DLOG_ERROR(logger)  LOG_NOTHING(LL_ERROR, logger)
#endif

// Note: Always logged.
#define DLOG_FATAL(logger)  LOG_FATAL(logger)

// Aliases - so we can do:
//      LOG(LL_FATAL, logger)
#define LOG_LL_TRACE(logger)    LOG_TRACE(logger)
#define LOG_LL_DEBUG(logger)    LOG_DEBUG(logger)
#define LOG_LL_INFO(logger)     LOG_INFO(logger)
#define LOG_LL_WARN(logger)     LOG_WARN(logger)
#define LOG_LL_ERROR(logger)    LOG_ERROR(logger)
#define LOG_LL_FATAL(logger)    LOG_FATAL(logger)

#define DLOG_LL_TRACE(logger)   DLOG_TRACE(logger)
#define DLOG_LL_DEBUG(logger)   DLOG_DEBUG(logger)
#define DLOG_LL_INFO(logger)    DLOG_INFO(logger)
#define DLOG_LL_WARN(logger)    DLOG_WARN(logger)
#define DLOG_LL_ERROR(logger)   DLOG_ERROR(logger)
#define DLOG_LL_FATAL(logger)   DLOG_FATAL(logger)


// Helper - if you want to do:
//      LOG(LL_FATAL, logger)
#define LOG(level, logger)  LOG_##level(logger)
#define DLOG(level, logger) DLOG_##level(logger)


// Log conditions.
#define LOG_IF(level, logger, condition)        !(condition) ? (void)0 : cpplog::helpers::VoidStreamClass() & LOG_##level(logger)
#define LOG_IF_NOT(level, logger, condition)    !!(condition) ? (void)0 : cpplog::helpers::VoidStreamClass() & LOG_##level(logger)

// Debug conditions.
#ifdef _DEBUG
#define DLOG_IF(level, logger, condition)       !(condition) ? (void)0 : cpplog::helpers::VoidStreamClass() & LOG_##level(logger)
#define DLOG_IF_NOT(level, logger, condition)   !!(condition) ? (void)0 : cpplog::helpers::VoidStreamClass() & LOG_##level(logger)
#else
#define DLOG_IF(level, logger, condition)       (true || !(condition)) ? (void)0 : \
                                                    cpplog::helpers::VoidStreamClass() & LOG_##level(logger)
#define DLOG_IF_NOT(level, logger, condition)       (true || !!(condition)) ? (void)0 : \
                                                    cpplog::helpers::VoidStreamClass() & LOG_##level(logger)
#endif


// Assertion helpers.
#define LOG_ASSERT(logger, condition)           LOG_IF_NOT(LL_FATAL, logger, (condition)) << "Assertion failed: " #condition
#define DLOG_ASSERT(logger, condition)          DLOG_IF_NOT(LL_FATAL, logger, (condition)) << "Assertion failed: " #condition


// Only include further helper macros if we are supposed to.
#ifdef CPPLOG_HELPER_MACROS

// The following CHECK_* functions act similar to a LOG_ASSERT, but with a bit more
// readability.

#define __CHECK(logger, condition, print)   LOG_IF(LL_FATAL, logger, !(condition)) << "Check failed: " print ": "

#define CHECK(logger, condition)            __CHECK(logger, condition, #condition)

#define CHECK_EQUAL(logger, ex1, ex2)       __CHECK(logger, (ex1) == (ex2), #ex1 " == " #ex2)
#define CHECK_LT(logger, ex1, ex2)          __CHECK(logger, (ex1) < (ex2), #ex1 " < " #ex2)
#define CHECK_GT(logger, ex1, ex2)          __CHECK(logger, (ex1) > (ex2), #ex1 " > " #ex2)
#define CHECK_LE(logger, ex1, ex2)          __CHECK(logger, (ex1) <= (ex2), #ex1 " <= " #ex2)
#define CHECK_GE(logger, ex1, ex2)          __CHECK(logger, (ex1) >= (ex2), #ex1 " >= " #ex2)

#define CHECK_NE(logger, ex1, ex2)          __CHECK(logger, (ex1) != (ex2), #ex1 " != " #ex2)
#define CHECK_NOT_EQUAL(logger, ex1, ex2)   __CHECK(logger, (ex1) != (ex2), #ex1 " != " #ex2)


// String helpers.
#define CHECK_STREQ(logger, s1, s2)         __CHECK(logger, strcmp((s1), (s2)) == 0, "") << s1 << " == " << s2
#define CHECK_STRNE(logger, s1, s2)         __CHECK(logger, strcmp((s1), (s2)) != 0, "") << s1 << " != " << s2

// NULL helpers.
#define CHECK_NULL(logger, exp)             __CHECK(logger, (exp) == NULL, #exp " == NULL")
#define CHECK_NOT_NULL(logger, exp)         __CHECK(logger, (exp) != NULL, #exp " != NULL")



// Debug versions of above.
#ifdef _DEBUG
#define DCHECK(logger, condition)           CHECK(logger, condition)
#define DCHECK_EQUAL(logger, ex1, ex2)      CHECK_EQUAL(logger, ex1, ex2)
#define DCHECK_LT(logger, ex1, ex2)         CHECK_LT(logger, ex1, ex2)
#define DCHECK_GT(logger, ex1, ex2)         CHECK_GT(logger, ex1, ex2)
#define DCHECK_LE(logger, ex1, ex2)         CHECK_LE(logger, ex1, ex2)
#define DCHECK_GE(logger, ex1, ex2)         CHECK_GE(logger, ex1, ex2)
#define DCHECK_NE(logger, ex1, ex2)         CHECK_NE(logger, ex1, ex2)
#define DCHECK_NOT_EQUAL(logger, ex1, ex2)  CHECK_NOT_EQUAL(logger, ex1, ex2)
#define DCHECK_STREQ(logger, s1, s2)        CHECK_STREQ(logger, s1, s2)
#define DCHECK_STRNE(logger, s1, s2)        CHECK_STRNE(logger, s1, s2)
#define DCHECK_NULL(logger, exp)            CHECK_NULL(logger, exp)
#define DCHECK_NOT_NULL(logger, exp)        CHECK_NOT_NULL(logger, exp)
#else
#define DCHECK(logger, condition)           while(false) CHECK(logger, condition)
#define DCHECK_EQUAL(logger, ex1, ex2)      while(false) CHECK_EQUAL(logger, ex1, ex2)
#define DCHECK_LT(logger, ex1, ex2)         while(false) CHECK_LT(logger, ex1, ex2)
#define DCHECK_GT(logger, ex1, ex2)         while(false) CHECK_GT(logger, ex1, ex2)
#define DCHECK_LE(logger, ex1, ex2)         while(false) CHECK_LE(logger, ex1, ex2)
#define DCHECK_GE(logger, ex1, ex2)         while(false) CHECK_GE(logger, ex1, ex2)
#define DCHECK_NE(logger, ex1, ex2)         while(false) CHECK_NE(logger, ex1, ex2)
#define DCHECK_NOT_EQUAL(logger, ex1, ex2)  while(false) CHECK_NOT_EQUAL(logger, ex1, ex2)
#define DCHECK_STREQ(logger, s1, s2)        while(false) CHECK_STREQ(logger, s1, s2)
#define DCHECK_STRNE(logger, s1, s2)        while(false) CHECK_STRNE(logger, s1, s2)
#define DCHECK_NULL(logger, exp)            while(false) CHECK_NULL(logger, exp)
#define DCHECK_NOT_NULL(logger, exp)        while(false) CHECK_NOT_NULL(logger, exp)
#endif


#endif

#endif //_CPPLOG_H
