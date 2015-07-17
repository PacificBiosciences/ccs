#pragma once
#ifndef _OUTPUT_DEBUG_STREAM_H
#define _OUTPUT_DEBUG_STREAM_H

#define NOMINMAX
#include <windows.h>
#undef DELETE

template <class _Elem>
class outputdebug_buf: public std::basic_stringbuf<_Elem>
{
private:
    inline void output_debug_string(const _Elem* e);

public:
    virtual int sync ( )
    {
        //outputdebug_buf<_Elem>::output_debug_string(std::basic_stringbuf<_Elem>::str().c_str());
        this->output_debug_string(std::basic_stringbuf<_Elem>::str().c_str());
        return 0;
    }
};

template<>
inline void outputdebug_buf<char>::output_debug_string(const char* e)
{
    ::OutputDebugStringA(e);
}

template<>
inline void outputdebug_buf<wchar_t>::output_debug_string(const wchar_t* e)
{
    ::OutputDebugStringW(e);
}

template <class _Elem>
class outputdebug_stream: public std::basic_ostream<_Elem>
{
public:
    outputdebug_stream() : std::basic_ostream<_Elem>(new outputdebug_buf<_Elem>())
    {}

    ~outputdebug_stream(){ delete this->rdbuf(); }
};

typedef outputdebug_stream<char> dbgwin_stream;
typedef outputdebug_stream<wchar_t> wdbgwin_stream;

#endif
