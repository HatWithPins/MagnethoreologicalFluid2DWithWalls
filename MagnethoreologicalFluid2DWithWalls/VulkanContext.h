#include <vulkan/vulkan.hpp>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <array>
#include <cstring>
#include <iostream>

class VulkanContext {
public:
    // Singleton access
    static VulkanContext& instance();

    // Non-copyable
    VulkanContext(const VulkanContext&) = delete;
    VulkanContext& operator=(const VulkanContext&) = delete;

    // Accessors
    VkDevice device() const;
    VkQueue computeQueue() const;
    uint32_t computeQueueFamily() const;
    VkDescriptorSetLayout descriptorSetLayout() const;
    VkPipelineLayout pipelineLayout() const;
    VkPhysicalDevice physicalDevice() const;
    VkPipeline forcesPipeline() const;
    VkPipeline distancesPipeline() const;
    VkPipeline sumPipeline() const;
    VkPipeline validationPipeline() const;

private:
    VulkanContext();
    ~VulkanContext();

    void createInstance();
    void pickPhysicalDevice();
    void createDevice();
    void createDescriptorSetLayout();
    void createPipelineLayout();
    void createPipelines();

    VkShaderModule loadShader(const char* path);
    VkPipeline createComputePipeline(const char* shaderPath);

    VkInstance m_instance{ VK_NULL_HANDLE };
    VkPhysicalDevice m_physicalDevice{ VK_NULL_HANDLE };
    VkDevice m_device{ VK_NULL_HANDLE };
    VkQueue m_computeQueue{ VK_NULL_HANDLE };
    uint32_t m_computeQueueFamily{ 0 };

    VkDescriptorSetLayout m_descriptorSetLayout{ VK_NULL_HANDLE };
    VkPipelineLayout m_pipelineLayout{ VK_NULL_HANDLE };

    VkPipeline m_pipelineDistances;
    VkPipeline m_pipelineForces;
    VkPipeline m_pipelineSum;
    VkPipeline m_pipelineValidation;
};