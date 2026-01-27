#include "SimulationContext.h"

static uint32_t findMemoryType(
    VkPhysicalDevice physicalDevice,
    uint32_t typeFilter,
    VkMemoryPropertyFlags properties)
{
    VkPhysicalDeviceMemoryProperties memProps;
    vkGetPhysicalDeviceMemoryProperties(physicalDevice, &memProps);

    for (uint32_t i = 0; i < memProps.memoryTypeCount; i++) {
        if ((typeFilter & (1 << i)) &&
            (memProps.memoryTypes[i].propertyFlags & properties) == properties) {
            return i;
        }
    }

    throw std::runtime_error("Failed to find suitable memory type");
}

SimulationContext::SimulationContext(size_t particleCount) {
    m_particles = particleCount;
    VulkanContext& vk = VulkanContext::instance();
    m_device = vk.device();

    createCommandPool();
    allocateCommandBuffer();
    createFence();
    createBuffers();
    createDescriptorSet();
    recordCommands();
}

SimulationContext::~SimulationContext() {
    if (m_device == VK_NULL_HANDLE) return;

    vkDestroyFence(m_device, m_fence, nullptr);

    for (int i = 0; i < 31; i++){
      vkDestroyBuffer(m_device, m_buffers[i], nullptr);
      vkFreeMemory(m_device, m_deviceMemory[i], nullptr);
    }

    vkDestroyCommandPool(m_device, m_commandPool, nullptr);
}

void SimulationContext::createCommandPool() {
    VkCommandPoolCreateInfo info{ VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO };
    info.queueFamilyIndex = VulkanContext::instance().computeQueueFamily();
    info.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;

    if (vkCreateCommandPool(m_device, &info, nullptr, &m_commandPool) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create command pool");
    }
}

void SimulationContext::allocateCommandBuffer() {
    VkCommandBufferAllocateInfo info{ VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO };
    info.commandPool = m_commandPool;
    info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    info.commandBufferCount = 4;

    if (vkAllocateCommandBuffers(m_device, &info, m_commandBuffers) != VK_SUCCESS) {
        throw std::runtime_error("Failed to allocate command buffer");
    }
}
void SimulationContext::recordCommand(int command, VkCommandBufferBeginInfo info, VkPipeline pipeline, VkPipelineLayout pipelineLayout, int numThreads){
    vkBeginCommandBuffer(m_commandBuffers[command], &info);

    vkCmdBindPipeline(m_commandBuffers[command], VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
    vkCmdBindDescriptorSets(m_commandBuffers[command], VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, 0, 1, &m_descriptorSet, 0, nullptr);

    VkMemoryBarrier hostToDevice{};
    hostToDevice.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
    hostToDevice.srcAccessMask = VK_ACCESS_HOST_WRITE_BIT;
    hostToDevice.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;

    vkCmdPipelineBarrier(
        m_commandBuffers[command],
        VK_PIPELINE_STAGE_HOST_BIT,
        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
        0,
        1, &hostToDevice,
        0, nullptr,
        0, nullptr
    );

    vkCmdDispatch(m_commandBuffers[command], numThreads, 1, 1);

    VkMemoryBarrier barrier{};
    barrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
    barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_HOST_READ_BIT;

    vkCmdPipelineBarrier(
        m_commandBuffers[command],
        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
        VK_PIPELINE_STAGE_HOST_BIT,
        0,
        1, &barrier,
        0, nullptr,
        0, nullptr
    );
    vkEndCommandBuffer(m_commandBuffers[command]);
}
void SimulationContext::recordCommands() {
    int matrixSize = m_particles * (m_particles - 1) / 2;;
    VulkanContext& vk = VulkanContext::instance();

    VkCommandBufferBeginInfo info{};
    info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    recordCommand(0, info, vk.forcesPipeline(),     vk.pipelineLayout(), matrixSize);
    recordCommand(1, info, vk.sumPipeline(),        vk.pipelineLayout(), m_particles);
    recordCommand(2, info, vk.distancesPipeline(),  vk.pipelineLayout(), matrixSize);
    recordCommand(3, info, vk.validationPipeline(), vk.pipelineLayout(), m_particles);
    
}

void SimulationContext::createFence() {
    VkFenceCreateInfo info{ VK_STRUCTURE_TYPE_FENCE_CREATE_INFO };
    info.flags = VK_FENCE_CREATE_SIGNALED_BIT;

    if (vkCreateFence(m_device, &info, nullptr, &m_fence) != VK_SUCCESS) {
        throw std::runtime_error("Failed to create fence");
    }
}

void SimulationContext::createBuffers() {
    /*
    31 buffers:
    0-2: x_0, y_0 and z_0
    3-5: x_1, y_1 and z_1
    6-8: F_x, F_y and F_z
    9: magnetic_field
    10: mode
    11: phase
    12: dimensions
    13: lenght
    14: field_direction
    15: particles
    16-18: delta_t, delta_t_original and t
    19-20: initial_idx_sum, last_idx_sum
    21-22: particle_0, particle_1
    23: matrix_size
    24: mason
    25: AR
    26: valid
    27-28: stress, stress_array
    29: r_array
    30: wall_velocity
    */
    int matrix_size = m_particles * (m_particles - 1) / 2;
    std::array<VkDeviceSize, 31> sizes = { 
        m_particles * sizeof(double), //x_0 
        m_particles * sizeof(double), //y_0
        m_particles * sizeof(double), //z_0
        m_particles * sizeof(double), //x_1
        m_particles * sizeof(double), //y_1
        m_particles * sizeof(double), //z_1
        matrix_size * sizeof(double), //F_x
        matrix_size * sizeof(double), //F_y
        matrix_size * sizeof(double), //F_z
        3 * sizeof(double),           //magnetic_field
        sizeof(int),                  //mode
        sizeof(int),                  //phase
        sizeof(int),                  //dimensions
        sizeof(int),                  //lenght
        sizeof(double),               //field_direction
        sizeof(int),                  //particles
        sizeof(double),               //delta_t
        sizeof(double),               //delta_t_original
        sizeof(double),               //t
        matrix_size * sizeof(int),    //initial_idx_sum
        matrix_size * sizeof(int),    //last_idx_sum
        matrix_size * sizeof(int),    //particle_0
        matrix_size * sizeof(int),    //particle_1
        sizeof(int),                  //matrix_size
        sizeof(double),               //mason
        sizeof(double),               //AR
        sizeof(int),                  //valid
        m_particles * sizeof(double), //stress_array
        sizeof(double),               //stress
        matrix_size * sizeof(double), //r_array
        sizeof(double)               //wall_velocity
    };
    for (int i = 0; i < 31; i++) {
        auto createBuffer = [&](VkBuffer& buffer, VkDeviceMemory& memory) {
            VkBufferCreateInfo bufInfo{ VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO };
            bufInfo.size = sizes[i];
            bufInfo.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;
            bufInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
            bufInfo.pQueueFamilyIndices = 0;

            if (vkCreateBuffer(m_device, &bufInfo, nullptr, &buffer) != VK_SUCCESS)
                throw std::runtime_error("Failed to create buffer");

            VkMemoryRequirements req;
            vkGetBufferMemoryRequirements(m_device, buffer, &req);

            VkMemoryAllocateInfo allocInfo{ VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO };
            allocInfo.allocationSize = req.size;
            allocInfo.memoryTypeIndex = findMemoryType(
                VulkanContext::instance().physicalDevice(),
                req.memoryTypeBits,
                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT
            );

            if (vkAllocateMemory(m_device, &allocInfo, nullptr, &memory) != VK_SUCCESS)
                throw std::runtime_error("Failed to allocate buffer memory");


        vkBindBufferMemory(m_device, buffer, memory, 0);

        };

    
        createBuffer(m_buffers[i], m_deviceMemory[i]);
    }
}

void SimulationContext::createDescriptorSet() {
    // Descriptor pool (one per simulation for simplicity)
    VkDescriptorPoolSize poolSize{};
    poolSize.type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    poolSize.descriptorCount = 31;

    VkDescriptorPoolCreateInfo poolInfo{ VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO };
    poolInfo.poolSizeCount = 1;
    poolInfo.pPoolSizes = &poolSize;
    poolInfo.maxSets = 1;

    VkDescriptorPool pool;
    vkCreateDescriptorPool(m_device, &poolInfo, nullptr, &pool);

    VkDescriptorSetAllocateInfo allocInfo{ VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO };
    allocInfo.descriptorPool = pool;
    allocInfo.descriptorSetCount = 1;

    VkDescriptorSetLayout layout =
        VulkanContext::instance().descriptorSetLayout();
    allocInfo.pSetLayouts = &layout;

    if (vkAllocateDescriptorSets(m_device, &allocInfo, &m_descriptorSet) != VK_SUCCESS)
        throw std::runtime_error("Failed to allocate descriptor set");

    std::array<VkWriteDescriptorSet, 31> writes{};
    std::array<VkDescriptorBufferInfo, 31> bufferInfos;

    for (int i = 0; i < 31; i++) {
        bufferInfos[i].buffer = m_buffers[i];
        bufferInfos[i].offset = 0;
        bufferInfos[i].range = VK_WHOLE_SIZE;

        writes[i] = { VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET };
        writes[i].dstSet = m_descriptorSet;
        writes[i].dstBinding = i;
        writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[i].descriptorCount = 1;
        writes[i].pBufferInfo = &bufferInfos[i];
    }

    vkUpdateDescriptorSets(m_device, static_cast<uint32_t>(writes.size()), writes.data(), 0, nullptr);
}

void SimulationContext::begin() {
    vkWaitForFences(m_device, 1, &m_fence, VK_TRUE, UINT64_MAX);
    vkResetFences(m_device, 1, &m_fence);
}

void SimulationContext::submit() {
    VulkanContext& vk = VulkanContext::instance();

    VkSubmitInfo submit{ VK_STRUCTURE_TYPE_SUBMIT_INFO };
    submit.commandBufferCount = 1;

    for (int i = 0; i < 4; i++) {
        submit.pCommandBuffers = &m_commandBuffers[i];

        vkQueueSubmit(
            vk.computeQueue(),
            1, &submit,
            m_fence
        );
        wait();
    }
}

void SimulationContext::wait() {
    int fenceWaitStatus = vkWaitForFences(m_device, 1, &m_fence, VK_TRUE, UINT64_MAX);
}

VkCommandBuffer SimulationContext::commandBuffer() const {
    return m_commandBuffers[0];
}

std::array<VkBuffer, 31> SimulationContext::hostBuffers() {
    return m_buffers;
}

VkDescriptorSet SimulationContext::descriptorSet() const {
    return m_descriptorSet;
}

void SimulationContext::SetX(double* x) {
    void* mapped = nullptr;
    vkMapMemory(m_device, m_deviceMemory[0], 0, m_particles * sizeof(double), 0, &mapped);

    double* xDevice = static_cast<double*>(mapped);
    for (int32_t I = 0; I < m_particles; ++I)
    {
        xDevice[I] = x[I];
    }
    vkUnmapMemory(m_device, m_deviceMemory[0]);
}
std::vector<double> SimulationContext::ReturnX() {
    std::vector<double> x(m_particles);
    void* mapped = nullptr;
    vkMapMemory(m_device, m_deviceMemory[0], 0, m_particles * sizeof(double), 0, &mapped);

    double* xDevice = static_cast<double*>(mapped);
    for (uint32_t i = 0; i < m_particles; ++i) {
        x[i] = xDevice[i];
    }

    vkUnmapMemory(m_device, m_deviceMemory[0]);

    return x;
}
void SimulationContext::SetY(double* y) {
    void* mapped = nullptr;
    vkMapMemory(m_device, m_deviceMemory[1], 0, m_particles * sizeof(double), 0, &mapped);

    double* yDevice = static_cast<double*>(mapped);
    for (int32_t I = 0; I < m_particles; ++I)
    {
        yDevice[I] = y[I];
    }
    vkUnmapMemory(m_device, m_deviceMemory[1]);
}
std::vector<double> SimulationContext::ReturnY() {
    std::vector<double> y(m_particles);
    void* mapped = nullptr;
    vkMapMemory(m_device, m_deviceMemory[1], 0, m_particles * sizeof(double), 0, &mapped);

    double* yDevice = static_cast<double*>(mapped);
    for (int32_t I = 0; I < m_particles; ++I)
    {
        y[I] = yDevice[I];
    }
    vkUnmapMemory(m_device, m_deviceMemory[1]);
    return y;
}
void SimulationContext::SetZ(double* z) {
    void* mapped = nullptr;
    vkMapMemory(m_device, m_deviceMemory[2], 0, m_particles * sizeof(double), 0, &mapped);

    double* zDevice = static_cast<double*>(mapped);
    for (int32_t I = 0; I < m_particles; ++I)
    {
        zDevice[I] = z[I];
    }
    vkUnmapMemory(m_device, m_deviceMemory[2]);
}
std::vector<double> SimulationContext::ReturnZ() {
    std::vector<double> z(m_particles);
    void* mapped = nullptr;
    vkMapMemory(m_device, m_deviceMemory[2], 0, m_particles * sizeof(double), 0, &mapped);

    double* zDevice = static_cast<double*>(mapped);
    for (int32_t I = 0; I < m_particles; ++I)
    {
        z[I] = zDevice[I];
    }
    vkUnmapMemory(m_device, m_deviceMemory[2]);
    return z;
}