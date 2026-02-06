#pragma once
#include <vulkan/vulkan.hpp>
#include <CL/cl.hpp>
#include <functional>
#include <condition_variable>

struct SubmitJob {
    VkSubmitInfo submit;
    VkFence fence;
};

struct SubmitJobOpenCL {
    std::function<void(cl::CommandQueue&)> record;
    std::shared_ptr<std::condition_variable> done;
    std::shared_ptr<std::mutex> mutex;
    std::shared_ptr<bool> finished;
};

struct StepSync {
    std::shared_ptr<std::mutex> mutex;
    std::shared_ptr<std::condition_variable> cv;
    std::shared_ptr<bool> finished;
};
