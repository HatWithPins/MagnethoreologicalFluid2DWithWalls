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

    void pushOpenCL(SubmitJobOpenCL job);
    bool popOpenCL(SubmitJobOpenCL& job);

private:
    std::queue<SubmitJob> queue;
    std::queue<SubmitJobOpenCL> queueOpenCL;
    std::mutex mutex;
    std::condition_variable cv;
    bool stop = false;
};