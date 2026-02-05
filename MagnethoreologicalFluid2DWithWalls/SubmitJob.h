#pragma once
#include <vulkan/vulkan.hpp>
#include <CL/cl.hpp>
#include <functional>

struct SubmitJob {
    VkSubmitInfo submit;
    VkFence fence;
};

struct SubmitJobOpenCL {
    std::function<void(cl::CommandQueue&)> record;
    std::condition_variable* done;
    bool* finished;
};
