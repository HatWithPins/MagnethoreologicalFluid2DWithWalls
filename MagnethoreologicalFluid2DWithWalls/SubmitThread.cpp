#include "SubmitThread.h"

extern SubmitQueue g_submitQueue;
extern std::atomic<bool> g_submitRunning;

void SubmitThreadFunc() {
    VulkanContext& vk = VulkanContext::instance();

    while (g_submitRunning) {
        SubmitJob job = g_submitQueue.pop();
        if (!g_submitRunning) break;

        vkQueueSubmit(
            vk.computeQueue(),
            1,
            &job.submit,
            job.fence
        );
    }
}

void SubmitThreadFuncOpenCL() {
    auto& queue = OpenCLContext::instance().queue();

    while (g_submitRunning) {
        SubmitJobOpenCL job;
        if (!g_submitQueue.popOpenCL(job)) {
            std::this_thread::sleep_for(std::chrono::microseconds(50));
            continue;
        }

        job.record(queue);
        queue.finish();

        {
            std::lock_guard<std::mutex> lock(*job.mutex);
            *job.finished = true;
        }
        job.done->notify_one();
    }
}
