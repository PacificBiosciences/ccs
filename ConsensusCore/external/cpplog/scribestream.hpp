#pragma once
#ifndef _SCRIBE_STREAM_H
#define _SCRIBE_STREAM_H

#include <protocol/TBinaryProtocol.h>
#include <transport/TSocket.h>
#include <transport/TTransportUtils.h>

#include "gen-cpp/scribe.h"

class scribe_buf : public std::basic_stringbuf<char>
{
private:
    scribe::thrift::scribeClient* m_client;
    boost::shared_ptr<apache::thrift::transport::TTransport> m_transport;

    std::string m_host;
    unsigned short m_port;
    std::string m_category;

public:
    scribe_buf() : m_client(NULL)
    {}

    virtual ~scribe_buf()
    {
        close();
    }

    bool open(std::string& host, unsigned short port, std::string& category, int timeout)
    {
        m_host = host;
        m_port = port;
        m_category = category;
        try
        {
            boost::shared_ptr<apache::thrift::transport::TSocket> socket(new apache::thrift::transport::TSocket(m_host, m_port));
            socket->setConnTimeout(timeout);
            socket->setRecvTimeout(timeout);
            socket->setSendTimeout(timeout);
            boost::shared_ptr<apache::thrift::transport::TTransport> transport(new apache::thrift::transport::TFramedTransport(socket));
            m_transport = transport;
            boost::shared_ptr<apache::thrift::protocol::TProtocol> proto(new apache::thrift::protocol::TBinaryProtocol(transport));
            m_client = new scribe::thrift::scribeClient(proto);
            transport->open();
        }
        catch (apache::thrift::TException& tx)
        {
            std::cerr << "Open scribe transport failed " << tx.what() << std::endl;
            return false;
        }
        return true;
    }

    virtual int sync()
    {
        if (m_client != NULL && m_transport && m_transport->isOpen())
        {
            scribe::thrift::LogEntry entry;
            entry.category = m_category.c_str();
            entry.message = str().c_str();
            std::vector<scribe::thrift::LogEntry> messages;
            messages.push_back(entry);
            try
            {
                int result = m_client->Log(messages);
                if (result != scribe::thrift::OK)
                {
                    std::cerr << "Log to scribe failed: " << result << " " << str().c_str() << std::endl;
                }
            }
            catch (apache::thrift::TException& e)
            {
                std::cerr << "Log to scribe exception: " << e.what() << " " << str().c_str() << std::endl;
            }
        }
        return 0;
    }

private:
    void close()
    {
        if (m_client != NULL)
        {
            try
            {
                m_transport->close();
            }
            catch (apache::thrift::TException& e) {}

            delete m_client;
            m_client = NULL;
        }
    }
};

class scribe_stream : public std::basic_ostream<char>
{
public:
    scribe_stream() : std::basic_ostream<char>(new scribe_buf()) { }
    ~scribe_stream() { delete rdbuf(); }

    void open(std::string host, unsigned short port, std::string category, int timeout)
    {
        scribe_buf * buf = (scribe_buf*)rdbuf();
        buf->open(host, port, category, timeout);
    }
};

#endif
