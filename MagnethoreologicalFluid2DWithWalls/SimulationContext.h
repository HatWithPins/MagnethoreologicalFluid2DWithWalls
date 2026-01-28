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
    SimulationContext(size_t particleCount, int dimensions, int length, double field_direction, double delta_t, double mason, double amplitude_relationship);
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
    void SetPhase(int phase);
    void SetMode(int mode);
    void SetMagneticField(double magnetic_field[3]);
    void SetMason(double mason);
    void SetWallVelocity(double wall_velocity);
    int ReturnValid();
    double ReturnTime();
    double ReturnDeltaT();
    double ReturnStress();

    VkCommandBuffer commandBuffer() const;
    VkDescriptorSet descriptorSet() const;

private:
    void createCommandPool();
    void allocateCommandBuffer();
    void createFence();
    void createBuffers();
    void preloadBuffers();
    void createDescriptorSet();
    void recordCommands();
    void recordCommand(int command, VkCommandBufferBeginInfo info, VkPipeline pipeline, VkPipelineLayout pipelineLayout, int numThreads);

    VkDevice m_device{ VK_NULL_HANDLE };

    VkCommandPool m_commandPool{ VK_NULL_HANDLE };
    VkCommandBuffer m_commandBuffers[4];
    VkFence m_fence{ VK_NULL_HANDLE };

    std::array<VkBuffer, 30> m_buffers;
    std::array<VkDeviceMemory, 30> m_deviceMemory;

    VkDescriptorSet m_descriptorSet{ VK_NULL_HANDLE };

    size_t m_particles;
    int m_dimensions;
    int m_length;
    double m_field_direction;
    double m_delta_t;
    double m_mason;
    double m_amplitude_relationship;
    int m_matrix_size;
    int m_valid;
};