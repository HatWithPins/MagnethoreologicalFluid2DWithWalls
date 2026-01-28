#pragma once
#include <vulkan/vulkan.hpp>

struct SubmitJob {
    VkSubmitInfo submit;
    VkFence fence;
};