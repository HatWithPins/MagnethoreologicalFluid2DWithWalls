#pragma once

#include "VulkanContext.h"
#include <vulkan/vulkan.hpp>
#include <vector>
#include <array>
#include <iostream>
#include <exception>
#include <stdexcept>

class SimulationContext {
public:
    SimulationContext(size_t particleCount);
    ~SimulationContext();

    // Non-copyable
    SimulationContext(const SimulationContext&) = delete;
    SimulationContext& operator=(const SimulationContext&) = delete;

    // Recording / execution
    void begin();
    void submit();
    void wait();
 
    void SetX(double* x);
    std::vector<double> ReturnX();
    void SetY(double* y);
    std::vector<double> ReturnY();
    void SetZ(double* z);
    std::vector<double> ReturnZ();

    VkCommandBuffer commandBuffer() const;

    // Buffers
    std::array<VkBuffer, 31> hostBuffers();

    VkDescriptorSet descriptorSet() const;

private:
    void createCommandPool();
    void allocateCommandBuffer();
    void createFence();
    void createBuffers();
    void createDescriptorSet();
    void recordCommands();
    void recordCommand(int command, VkCommandBufferBeginInfo info, VkPipeline pipeline, VkPipelineLayout pipelineLayout, int numThreads);

    VkDevice m_device{ VK_NULL_HANDLE };

    VkCommandPool m_commandPool{ VK_NULL_HANDLE };
    VkCommandBuffer m_commandBuffers[4];
    VkFence m_fence{ VK_NULL_HANDLE };

    std::array<VkBuffer, 31> m_buffers;
    std::array<VkDeviceMemory, 31> m_deviceMemory;

    VkDescriptorSet m_descriptorSet{ VK_NULL_HANDLE };

    size_t m_particles;
};