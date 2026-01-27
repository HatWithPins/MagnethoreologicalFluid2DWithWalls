#include "VulkanContext.h"

VulkanContext& VulkanContext::instance() {
    static VulkanContext ctx;
    return ctx;
}

VulkanContext::VulkanContext() {
    createInstance();
    pickPhysicalDevice();
    createDevice();
    createDescriptorSetLayout();
    createPipelineLayout();
    createPipelines();
}

VulkanContext::~VulkanContext() {
    if (m_pipelineLayout)
        vkDestroyPipelineLayout(m_device, m_pipelineLayout, nullptr);

    if (m_descriptorSetLayout)
        vkDestroyDescriptorSetLayout(m_device, m_descriptorSetLayout, nullptr);

    if (m_device)
        vkDestroyDevice(m_device, nullptr);

    if (m_instance)
        vkDestroyInstance(m_instance, nullptr);
}

void VulkanContext::createInstance() {
    VkApplicationInfo appInfo{ VK_STRUCTURE_TYPE_APPLICATION_INFO };
    appInfo.pApplicationName = "ComputeSimulation";
    appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.pEngineName = "None";
    appInfo.engineVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.apiVersion = VK_API_VERSION_1_2;

    VkInstanceCreateInfo info{ VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO };
    info.pApplicationInfo = &appInfo;

    if (vkCreateInstance(&info, nullptr, &m_instance) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create Vulkan instance");
    }
}

void VulkanContext::pickPhysicalDevice() {
    uint32_t count = 0;
    vkEnumeratePhysicalDevices(m_instance, &count, nullptr);
    if (count == 0)
        throw std::runtime_error("No Vulkan devices found");

    std::vector<VkPhysicalDevice> devices(count);
    vkEnumeratePhysicalDevices(m_instance, &count, devices.data());

    for (VkPhysicalDevice dev : devices) {
        uint32_t queueCount = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(dev, &queueCount, nullptr);

        std::vector<VkQueueFamilyProperties> props(queueCount);
        vkGetPhysicalDeviceQueueFamilyProperties(dev, &queueCount, props.data());

        for (uint32_t i = 0; i < queueCount; ++i) {
            if (props[i].queueFlags & VK_QUEUE_COMPUTE_BIT) {
                m_physicalDevice = dev;
                m_computeQueueFamily = i;
                return;
            }
        }
    }

    throw std::runtime_error("No compute-capable queue found");
}

void VulkanContext::createDevice() {
    float priority = 1.0f;

    VkDeviceQueueCreateInfo queueInfo{ VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO };
    queueInfo.queueFamilyIndex = m_computeQueueFamily;
    queueInfo.queueCount = 1;
    queueInfo.pQueuePriorities = &priority;

    VkPhysicalDeviceFeatures features{};
    features.shaderFloat64 = VK_TRUE;

    VkDeviceCreateInfo info{ VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO };
    info.queueCreateInfoCount = 1;
    info.pQueueCreateInfos = &queueInfo;
    info.pEnabledFeatures = &features;

    if (vkCreateDevice(m_physicalDevice, &info, nullptr, &m_device) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create Vulkan device");
    }

    vkGetDeviceQueue(m_device, m_computeQueueFamily, 0, &m_computeQueue);
}

void VulkanContext::createDescriptorSetLayout() {
    const int num_bindings = 31; //One binding for each original OpenCL buffer.
    std::array<VkDescriptorSetLayoutBinding, num_bindings> bindings{};

    for (int i = 0; i < num_bindings; i++) {
        bindings[i].binding = i;
        bindings[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[i].descriptorCount = 1;
        bindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    }

    VkDescriptorSetLayoutCreateInfo info{
        VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO
    };
    info.bindingCount = static_cast<uint32_t>(bindings.size());
    info.pBindings = bindings.data();

    if (vkCreateDescriptorSetLayout(m_device, &info, nullptr, &m_descriptorSetLayout) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create descriptor set layout");
    }
}

void VulkanContext::createPipelineLayout() {
    VkPipelineLayoutCreateInfo info{
        VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO
    };
    info.setLayoutCount = 1;
    info.pSetLayouts = &m_descriptorSetLayout;

    if (vkCreatePipelineLayout(m_device, &info, nullptr, &m_pipelineLayout) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create pipeline layout");
    }
}

VkShaderModule VulkanContext::loadShader(const char* path) {
    std::vector<char> ShaderContents;
    if (std::ifstream ShaderFile{path, std::ios::binary | std::ios::ate })
    {
        const size_t FileSize = ShaderFile.tellg();
        ShaderFile.seekg(0);
        ShaderContents.resize(FileSize, '\0');
        ShaderFile.read(ShaderContents.data(), FileSize);
    }

    VkShaderModuleCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    info.codeSize = ShaderContents.size();
    info.pCode = reinterpret_cast<const uint32_t*>(ShaderContents.data());

    VkShaderModule module;
    vkCreateShaderModule(m_device, &info, nullptr, &module);
    return module;
}

VkPipeline VulkanContext::createComputePipeline(const char* shaderPath) {
    VkShaderModule shader = loadShader(shaderPath);

    VkPipelineShaderStageCreateInfo stage{};
    stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    stage.module = shader;
    stage.pName = "Main"; //This assumes that all shaders have the same entry point.

    VkComputePipelineCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
    info.stage = stage;
    info.layout = m_pipelineLayout;

    VkPipeline pipeline;
    int pipelineCode = vkCreateComputePipelines(m_device, VK_NULL_HANDLE, 1, &info, nullptr, &pipeline);

    vkDestroyShaderModule(m_device, shader, nullptr);
    return pipeline;
}

void VulkanContext::createPipelines() {
    m_pipelineDistances  = createComputePipeline("../shaders/distances.spv");
    m_pipelineForces     = createComputePipeline("../shaders/forces.spv");
    m_pipelineSum        = createComputePipeline("../shaders/sum.spv");
    m_pipelineValidation = createComputePipeline("../shaders/validation.spv");
}

VkDevice VulkanContext::device() const {
    return m_device;
}

VkQueue VulkanContext::computeQueue() const {
    return m_computeQueue;
}

uint32_t VulkanContext::computeQueueFamily() const {
    return m_computeQueueFamily;
}

VkDescriptorSetLayout VulkanContext::descriptorSetLayout() const {
    return m_descriptorSetLayout;
}

VkPipelineLayout VulkanContext::pipelineLayout() const {
    return m_pipelineLayout;
}

VkPhysicalDevice VulkanContext::physicalDevice() const {
    return m_physicalDevice;
}
VkPipeline VulkanContext::forcesPipeline() const {
    return m_pipelineForces;
}
VkPipeline VulkanContext::distancesPipeline() const {
    return m_pipelineDistances;
}
VkPipeline VulkanContext::sumPipeline() const {
    return m_pipelineSum;
}
VkPipeline VulkanContext::validationPipeline() const {
    return m_pipelineValidation;
}