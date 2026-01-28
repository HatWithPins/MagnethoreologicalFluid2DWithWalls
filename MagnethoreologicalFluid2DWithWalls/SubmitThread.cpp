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