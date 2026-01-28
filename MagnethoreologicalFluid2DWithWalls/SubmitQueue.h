#pragma once
#include <queue>
#include <mutex>
#include <condition_variable>
#include "SubmitJob.h"

class SubmitQueue {
public:
    void push(SubmitJob job);
    SubmitJob pop();
    void shutdown();

private:
    std::queue<SubmitJob> queue;
    std::mutex mutex;
    std::condition_variable cv;
    bool stop = false;
};