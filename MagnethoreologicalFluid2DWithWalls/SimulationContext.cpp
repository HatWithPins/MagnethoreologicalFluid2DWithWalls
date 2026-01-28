#include "SimulationContext.h"

extern SubmitQueue g_submitQueue;

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

SimulationContext::SimulationContext(size_t particleCount, int dimensions, int length, double field_direction, double delta_t, double mason, double amplitude_relationship) {
    m_particles = particleCount;
    m_dimensions = dimensions;
    m_length = length;
    m_field_direction = field_direction;
    m_delta_t = delta_t;
    m_mason = mason;
    m_amplitude_relationship = amplitude_relationship;
    m_matrix_size = m_particles * (m_particles - 1) / 2;
    VulkanContext& vk = VulkanContext::instance();
    m_device = vk.device();

    createCommandPool();
    allocateCommandBuffer();
    createFence();
    createBuffers();
    createDescriptorSet();
    preloadBuffers();
    recordCommands();
}

SimulationContext::~SimulationContext() {
    if (m_device == VK_NULL_HANDLE) return;

    vkDestroyFence(m_device, m_fence, nullptr);

    for (int i = 0; i < 30; i++){
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
    info.commandBufferCount = 2;

    if (vkAllocateCommandBuffers(m_device, &info, m_commandBuffers) != VK_SUCCESS) {
        throw std::runtime_error("Failed to allocate command buffer");
    }
}
void SimulationContext::recordCommand(int command, VkCommandBufferBeginInfo info, VkPipeline pipeline, VkPipelineLayout pipelineLayout, int numThreads){
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
}
void SimulationContext::recordCommands() {
    int matrixSize = m_particles * (m_particles - 1) / 2;;
    VulkanContext& vk = VulkanContext::instance();

    int limit = 32;

    VkCommandBufferBeginInfo info{};
    info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

    vkBeginCommandBuffer(m_commandBuffers[0], &info);
    recordCommand(0, info, vk.validationPipeline(), vk.pipelineLayout(), (m_particles + limit - 1) / limit);
    recordCommand(0, info, vk.forcesPipeline(),     vk.pipelineLayout(), (matrixSize+limit-1) / limit);
    recordCommand(0, info, vk.sumPipeline(),        vk.pipelineLayout(), (m_particles+limit -1) / limit);
    recordCommand(0, info, vk.distancesPipeline(),  vk.pipelineLayout(), (matrixSize + limit -1) / limit);
    vkEndCommandBuffer(m_commandBuffers[0]);
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
    30 buffers:
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
    29: wall_velocity
    */
    int matrix_size = m_particles * (m_particles - 1) / 2;
    std::array<VkDeviceSize, 30> sizes = { 
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
        sizeof(double)               //wall_velocity
    };
    for (int i = 0; i < 30; i++) {
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
void SimulationContext::preloadBuffers(){
    int* initial_indices_sum = new int[m_particles];
    int* last_indices_sum = new int[m_particles];
    int* particle_0 = new int[m_matrix_size];
    int* particle_1 = new int[m_matrix_size];

    for (int i = 0; i < m_particles; i++) {
        initial_indices_sum[i] = 0;
        last_indices_sum[i] = 0;
    }
    for (int i = 1; i < m_particles; i++) {
        for (int j = 1; j <= i; j++) {
            initial_indices_sum[i] += m_particles - j;
        }
    }
    last_indices_sum[0] = m_particles - 1;
    for (int i = 1; i < m_particles; i++) {
        for (int j = 1; j <= i + 1; j++) {
            last_indices_sum[i] += m_particles - j;
        }
    }

    int index = 0;
    for (int i = 0; i < m_particles - 1; i++) {
        for (int j = i + 1; j < m_particles; j++) {
            particle_1[index] = i;
            particle_0[index] = j;
            index++;
        }
    }

    void* mapped = nullptr;
    vkMapMemory(m_device, m_deviceMemory[19], 0, m_particles * sizeof(int), 0, &mapped);
    int* bufferDevice = static_cast<int*>(mapped);
    for (int32_t I = 0; I < m_particles; ++I)
    {
        bufferDevice[I] = initial_indices_sum[I];
    }
    vkUnmapMemory(m_device, m_deviceMemory[19]);


    vkMapMemory(m_device, m_deviceMemory[20], 0, m_particles * sizeof(int), 0, &mapped);
    bufferDevice = static_cast<int*>(mapped);
    for (int32_t I = 0; I < m_particles; ++I)
    {
        bufferDevice[I] = last_indices_sum[I];
    }
    vkUnmapMemory(m_device, m_deviceMemory[20]);

    vkMapMemory(m_device, m_deviceMemory[21], 0, m_matrix_size * sizeof(int), 0, &mapped);
    bufferDevice = static_cast<int*>(mapped);
    for (int32_t I = 0; I < m_matrix_size; ++I)
    {
        bufferDevice[I] = particle_0[I];
    }
    vkUnmapMemory(m_device, m_deviceMemory[21]);

    vkMapMemory(m_device, m_deviceMemory[22], 0, m_matrix_size * sizeof(int), 0, &mapped);
    bufferDevice = static_cast<int*>(mapped);
    for (int32_t I = 0; I < m_matrix_size; ++I)
    {
        bufferDevice[I] = particle_1[I];
    }
    vkUnmapMemory(m_device, m_deviceMemory[22]);

    vkMapMemory(m_device, m_deviceMemory[9], 0, 3 * sizeof(double), 0, &mapped);
    double* bufferDeviceDouble = static_cast<double*>(mapped);
    bufferDeviceDouble[0] = 0.0;
    bufferDeviceDouble[1] = (m_dimensions == 2);
    bufferDeviceDouble[2] = (m_dimensions == 3);
    vkUnmapMemory(m_device, m_deviceMemory[9]);

    vkMapMemory(m_device, m_deviceMemory[12], 0,  sizeof(int), 0, &mapped);
    bufferDevice = static_cast<int*>(mapped);
    bufferDevice[0] = m_dimensions;
    vkUnmapMemory(m_device, m_deviceMemory[12]);

    vkMapMemory(m_device, m_deviceMemory[13], 0, sizeof(int), 0, &mapped);
    bufferDevice = static_cast<int*>(mapped);
    bufferDevice[0] = m_length;
    vkUnmapMemory(m_device, m_deviceMemory[13]);

    vkMapMemory(m_device, m_deviceMemory[14], 0, sizeof(double), 0, &mapped);
    bufferDeviceDouble = static_cast<double*>(mapped);
    bufferDeviceDouble[0] = m_field_direction;
    vkUnmapMemory(m_device, m_deviceMemory[14]);

    vkMapMemory(m_device, m_deviceMemory[15], 0, sizeof(int), 0, &mapped);
    bufferDevice = static_cast<int*>(mapped);
    bufferDevice[0] = m_particles;
    vkUnmapMemory(m_device, m_deviceMemory[15]);

    vkMapMemory(m_device, m_deviceMemory[16], 0, sizeof(double), 0, &mapped);
    bufferDeviceDouble = static_cast<double*>(mapped);
    bufferDeviceDouble[0] = m_delta_t;
    vkUnmapMemory(m_device, m_deviceMemory[16]);

    vkMapMemory(m_device, m_deviceMemory[17], 0, sizeof(double), 0, &mapped);
    bufferDeviceDouble = static_cast<double*>(mapped);
    bufferDeviceDouble[0] = m_delta_t;
    vkUnmapMemory(m_device, m_deviceMemory[17]);

    vkMapMemory(m_device, m_deviceMemory[18], 0, sizeof(double), 0, &mapped);
    bufferDeviceDouble = static_cast<double*>(mapped);
    bufferDeviceDouble[0] = 0.0;
    vkUnmapMemory(m_device, m_deviceMemory[18]);

    vkMapMemory(m_device, m_deviceMemory[23], 0, sizeof(int), 0, &mapped);
    bufferDevice = static_cast<int*>(mapped);
    bufferDevice[0] = m_matrix_size;
    vkUnmapMemory(m_device, m_deviceMemory[23]);

    vkMapMemory(m_device, m_deviceMemory[24], 0, sizeof(double), 0, &mapped);
    bufferDeviceDouble = static_cast<double*>(mapped);
    bufferDeviceDouble[0] = m_mason;
    vkUnmapMemory(m_device, m_deviceMemory[24]);

    vkMapMemory(m_device, m_deviceMemory[25], 0, sizeof(double), 0, &mapped);
    bufferDeviceDouble = static_cast<double*>(mapped);
    bufferDeviceDouble[0] = m_amplitude_relationship;
    vkUnmapMemory(m_device, m_deviceMemory[25]);
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

    std::array<VkWriteDescriptorSet, 30> writes{};
    std::array<VkDescriptorBufferInfo, 30> bufferInfos;

    for (int i = 0; i < 30; i++) {
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
    SubmitJob job{};
    job.submit.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    job.submit.commandBufferCount = 1;
    job.submit.pCommandBuffers = &m_commandBuffers[0];
    job.fence = m_fence;

    g_submitQueue.push(job);
    
    wait();

    void* mapped = nullptr;
    vkMapMemory(m_device, m_deviceMemory[26], 0, sizeof(int), 0, &mapped);
    int* bufferDevice = static_cast<int*>(mapped);
    m_valid = bufferDevice[0];
    vkUnmapMemory(m_device, m_deviceMemory[26]);
}

void SimulationContext::wait() {
    int fenceWaitStatus = vkWaitForFences(m_device, 1, &m_fence, VK_TRUE, UINT64_MAX);
}

VkCommandBuffer SimulationContext::commandBuffer() const {
    return m_commandBuffers[0];
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
void SimulationContext::SetPhase(int phase){
    void* mapped = nullptr;

    vkMapMemory(m_device, m_deviceMemory[11], 0, sizeof(int), 0, &mapped);
    int* bufferDevice = static_cast<int*>(mapped);
    bufferDevice[0] = phase;
    vkUnmapMemory(m_device, m_deviceMemory[11]);
}
void SimulationContext::SetMode(int mode){
    void* mapped = nullptr;

    vkMapMemory(m_device, m_deviceMemory[10], 0, sizeof(int), 0, &mapped);
    int* bufferDevice = static_cast<int*>(mapped);
    bufferDevice[0] = mode;
    vkUnmapMemory(m_device, m_deviceMemory[10]);
}
void SimulationContext::SetMagneticField(double magnetic_field[3]){
    void* mapped = nullptr;

    vkMapMemory(m_device, m_deviceMemory[9], 0, 3 * sizeof(double), 0, &mapped);
    double* bufferDevice = static_cast<double*>(mapped);
    bufferDevice[0] = magnetic_field[0];
    bufferDevice[1] = magnetic_field[1];
    bufferDevice[2] = magnetic_field[2];
    vkUnmapMemory(m_device, m_deviceMemory[9]);
}
void SimulationContext::SetMason(double mason){
    void* mapped = nullptr;

    vkMapMemory(m_device, m_deviceMemory[24], 0, sizeof(double), 0, &mapped);
    double* bufferDevice = static_cast<double*>(mapped);
    bufferDevice[0] = mason;
    vkUnmapMemory(m_device, m_deviceMemory[24]);
}
void SimulationContext::SetWallVelocity(double wall_velocity){
    void* mapped = nullptr;

    vkMapMemory(m_device, m_deviceMemory[29], 0, sizeof(double), 0, &mapped);
    double* bufferDevice = static_cast<double*>(mapped);
    bufferDevice[0] = wall_velocity;
    vkUnmapMemory(m_device, m_deviceMemory[29]);
}
int SimulationContext::ReturnValid(){
    return m_valid;
}
double SimulationContext::ReturnTime(){
    void* mapped = nullptr;

    vkMapMemory(m_device, m_deviceMemory[18], 0, sizeof(double), 0, &mapped);
    double* bufferDevice = static_cast<double*>(mapped);
    double time = bufferDevice[0];
    vkUnmapMemory(m_device, m_deviceMemory[18]);

    return time;
}
double SimulationContext::ReturnDeltaT(){
    void* mapped = nullptr;

    vkMapMemory(m_device, m_deviceMemory[16], 0, sizeof(double), 0, &mapped);
    double* bufferDevice = static_cast<double*>(mapped);
    double delta_t = bufferDevice[0];
    vkUnmapMemory(m_device, m_deviceMemory[16]);

    return delta_t;
}
double SimulationContext::ReturnStress(){
    void* mapped = nullptr;

    vkMapMemory(m_device, m_deviceMemory[27], 0, sizeof(double), 0, &mapped);
    double* bufferDevice = static_cast<double*>(mapped);
    double stress = bufferDevice[0];
    vkUnmapMemory(m_device, m_deviceMemory[27]);

    return stress;
}