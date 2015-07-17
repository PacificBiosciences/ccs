#pragma once
#ifndef _CONCURRENT_QUEUE_H
#define _CONCURRENT_QUEUE_H

#include <queue>
#include <boost/thread.hpp>

template<typename Data>
class concurrent_queue
{
private:
    std::queue<Data> the_queue;
    mutable boost::mutex the_mutex;
    boost::condition_variable the_condition_variable;
public:
    void push(Data const& data)
    {
        boost::lock_guard<boost::mutex> lock(the_mutex);
        the_queue.push(data);
        the_condition_variable.notify_one();
    }

    bool empty() const
    {
        boost::lock_guard<boost::mutex> lock(the_mutex);
        return the_queue.empty();
    }

    bool try_pop(Data& popped_value)
    {
        boost::unique_lock<boost::mutex> lock(the_mutex);
        if( the_queue.empty() )
        {
            return false;
        }

        popped_value = the_queue.front();
        the_queue.pop();
        return true;
    }

    void wait_and_pop(Data& popped_value)
    {
        boost::unique_lock<boost::mutex> lock(the_mutex);

        while( the_queue.empty() )
        {
            the_condition_variable.wait(lock);
        }

        popped_value = the_queue.front();
        the_queue.pop();
    }

};

#endif
