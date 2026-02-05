#pragma once
#include <CL/cl.hpp>

class OpenCLContext {
public:
    static OpenCLContext& instance();

    cl::Context& context();
    cl::Device& device();
    cl::CommandQueue& queue();

private:
    OpenCLContext();
    cl::Platform m_platform;
    cl::Device   m_device;
    cl::Context  m_context;
    cl::CommandQueue m_queue;
};
