#include "OpenCLContext.h"

OpenCLContext& OpenCLContext::instance() {
    static OpenCLContext ctx;
    return ctx;
}

OpenCLContext::OpenCLContext() {
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    m_platform = platforms[0];

    std::vector<cl::Device> devices;
    m_platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
    m_device = devices[0];

    m_context = cl::Context({ m_device });
    m_queue = cl::CommandQueue(m_context, m_device);
}

cl::Context& OpenCLContext::context() { return m_context; }
cl::Device& OpenCLContext::device() { return m_device; }
cl::CommandQueue& OpenCLContext::queue() { return m_queue; }
