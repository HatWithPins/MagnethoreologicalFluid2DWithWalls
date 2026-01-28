#include "SubmitQueue.h"

void SubmitQueue::push(SubmitJob job) {
    {
        std::lock_guard<std::mutex> lock(mutex);
        queue.push(std::move(job));
    }
    cv.notify_one();
}

SubmitJob SubmitQueue::pop() {
    std::unique_lock<std::mutex> lock(mutex);
    cv.wait(lock, [&] { return !queue.empty() || stop; });

    SubmitJob job = std::move(queue.front());
    queue.pop();
    return job;
}

void SubmitQueue::shutdown() {
    stop = true;
    cv.notify_all();
}