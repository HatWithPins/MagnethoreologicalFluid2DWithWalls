#include "SimulationContextOCL.h"
#include "OpenCLContext.h"
#include <fstream>
#include <vector>
#include <stdexcept>

static cl::Program loadProgram(
    const char* path,
    cl::Context& ctx,
    cl::Device& dev)
{
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Failed to open kernel file");

    std::string src(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>()
    );

    cl::Program prog(ctx, src);
    prog.build({ dev });
    return prog;
}

SimulationContextOCL::SimulationContextOCL(int particles, int dimensions, int length, double field_direction, double ma, double ar, double original_delta_t)
    : m_particles(particles),
    m_dimensions(dimensions),
    m_length(length),
    m_field_direction(field_direction),
    m_ar(ar),
    m_mason(ma),
    m_delta_t(original_delta_t),
    m_context(OpenCLContext::instance().context()),
    m_device(OpenCLContext::instance().device()),
    m_queue(OpenCLContext::instance().queue())
{
    createBuffers();
    loadKernels();
}

void SimulationContextOCL::createBuffers() {
    size_t pbytes = m_particles * sizeof(double);
    size_t fcount = m_particles * (m_particles - 1) / 2;
    m_matrix_size = fcount;
    double time = 0;

    buf_x0 = cl::Buffer(m_context, CL_MEM_READ_WRITE, pbytes);
    buf_y0 = cl::Buffer(m_context, CL_MEM_READ_WRITE, pbytes);
    buf_z0 = cl::Buffer(m_context, CL_MEM_READ_WRITE, pbytes);

    buf_x1 = cl::Buffer(m_context, CL_MEM_READ_WRITE, pbytes);
    buf_y1 = cl::Buffer(m_context, CL_MEM_READ_WRITE, pbytes);
    buf_z1 = cl::Buffer(m_context, CL_MEM_READ_WRITE, pbytes);

    buf_fx = cl::Buffer(m_context, CL_MEM_READ_WRITE, fcount * sizeof(double));
    buf_fy = cl::Buffer(m_context, CL_MEM_READ_WRITE, fcount * sizeof(double));
    buf_fz = cl::Buffer(m_context, CL_MEM_READ_WRITE, fcount * sizeof(double));

    buf_ar               = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(double));
    buf_mason            = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(double));
    buf_magnetic_field   = cl::Buffer(m_context, CL_MEM_READ_WRITE, 3 * sizeof(double));
    buf_field_direction  = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(double));
    buf_wall_velocity    = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(double));
    buf_mode             = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(int));;
    buf_length           = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(int));
    buf_phase            = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(int));
    buf_dimensions       = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(int));
    buf_particles        = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(int));
    buf_matrix_size      = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(int));
    buf_t                = cl::Buffer(m_context, CL_MEM_READ_WRITE, sizeof(double));
    buf_delta_t          = cl::Buffer(m_context, CL_MEM_READ_WRITE, sizeof(double));
    buf_original_delta_t = cl::Buffer(m_context, CL_MEM_READ_ONLY, sizeof(double));
    buf_valid            = cl::Buffer(m_context, CL_MEM_READ_WRITE, sizeof(int));
    buf_initial_idx      = cl::Buffer(m_context, CL_MEM_READ_ONLY, fcount * sizeof(int));
    buf_last_idx         = cl::Buffer(m_context, CL_MEM_READ_ONLY, fcount * sizeof(int));
    buf_particle_0       = cl::Buffer(m_context, CL_MEM_READ_ONLY, fcount * sizeof(int));
    buf_particle_1       = cl::Buffer(m_context, CL_MEM_READ_ONLY, fcount * sizeof(int));
    buf_stress           = cl::Buffer(m_context, CL_MEM_READ_WRITE, sizeof(double));
    buf_stress_array     = cl::Buffer(m_context, CL_MEM_READ_WRITE, pbytes);

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

    double magnetic_field[3];
    if (m_dimensions == 3) {
        magnetic_field[0] = 0.0;
        magnetic_field[1] = 0.0;
        magnetic_field[2] = 1.0;
    }
    else {
        magnetic_field[0] = 0.0;
        magnetic_field[1] = 1.0;
        magnetic_field[2] = 0.0;
    }

    m_queue.enqueueWriteBuffer(buf_particles, CL_TRUE, 0, sizeof(int), &m_particles);
    m_queue.enqueueWriteBuffer(buf_dimensions, CL_TRUE, 0, sizeof(int), &m_dimensions);
    m_queue.enqueueWriteBuffer(buf_length, CL_TRUE, 0, sizeof(int), &m_length);
    m_queue.enqueueWriteBuffer(buf_field_direction, CL_TRUE, 0, sizeof(double), &m_field_direction);
    m_queue.enqueueWriteBuffer(buf_matrix_size, CL_TRUE, 0, sizeof(int), &m_matrix_size);
    m_queue.enqueueWriteBuffer(buf_delta_t, CL_TRUE, 0, sizeof(double), &m_delta_t);
    m_queue.enqueueWriteBuffer(buf_original_delta_t, CL_TRUE, 0, sizeof(double), &m_delta_t);
    m_queue.enqueueWriteBuffer(buf_mason, CL_TRUE, 0, sizeof(double), &m_mason);
    m_queue.enqueueWriteBuffer(buf_ar, CL_TRUE, 0, sizeof(double), &m_ar);
    m_queue.enqueueWriteBuffer(buf_t, CL_TRUE, 0, sizeof(double), &time);
    m_queue.enqueueWriteBuffer(buf_magnetic_field, CL_TRUE, 0, 3 * sizeof(double), &magnetic_field);
    m_queue.enqueueWriteBuffer(buf_particle_0, CL_TRUE, 0, m_matrix_size * sizeof(int), particle_0);
    m_queue.enqueueWriteBuffer(buf_particle_1, CL_TRUE, 0, m_matrix_size * sizeof(int), particle_1);
    m_queue.enqueueWriteBuffer(buf_initial_idx, CL_TRUE, 0, m_particles * sizeof(int), initial_indices_sum);
    m_queue.enqueueWriteBuffer(buf_last_idx, CL_TRUE, 0, m_particles * sizeof(int), last_indices_sum);
}

void SimulationContextOCL::loadKernels() {
    prog_forces    = loadProgram("forces.cl", m_context, m_device);
    prog_sum       = loadProgram("sum.cl", m_context, m_device);
    prog_validate  = loadProgram("validation.cl", m_context, m_device);
    prog_distances = loadProgram("distances.cl", m_context, m_device);

    kernel_forces    = cl::Kernel(prog_forces, "forces");
    kernel_sum       = cl::Kernel(prog_sum, "integrate");
    kernel_validate  = cl::Kernel(prog_validate, "validation");
    kernel_distances = cl::Kernel(prog_distances, "distances");

    kernel_forces.setArg(0, buf_x0);
    kernel_forces.setArg(1, buf_y0);
    kernel_forces.setArg(2, buf_z0);
    kernel_forces.setArg(3, buf_dimensions);
    kernel_forces.setArg(4, buf_magnetic_field);
    kernel_forces.setArg(5, buf_length);
    kernel_forces.setArg(6, buf_particles);
    kernel_forces.setArg(7, buf_matrix_size);
    kernel_forces.setArg(8, buf_particle_0);
    kernel_forces.setArg(9, buf_particle_1);
    kernel_forces.setArg(10, buf_fx);
    kernel_forces.setArg(11, buf_fy);
    kernel_forces.setArg(12, buf_fz);

    kernel_sum.setArg(0, buf_x0);
    kernel_sum.setArg(1, buf_y0);
    kernel_sum.setArg(2, buf_z0);
    kernel_sum.setArg(3, buf_dimensions);
    kernel_sum.setArg(4, buf_length);
    kernel_sum.setArg(5, buf_particles);
    kernel_sum.setArg(6, buf_delta_t);
    kernel_sum.setArg(7, buf_fx);
    kernel_sum.setArg(8, buf_fy);
    kernel_sum.setArg(9, buf_fz);
    kernel_sum.setArg(10, buf_initial_idx);
    kernel_sum.setArg(11, buf_last_idx);
    kernel_sum.setArg(12, buf_matrix_size);
    kernel_sum.setArg(13, buf_valid);
    kernel_sum.setArg(14, buf_x1);
    kernel_sum.setArg(15, buf_y1);
    kernel_sum.setArg(16, buf_z1);
    kernel_sum.setArg(17, buf_mode);
    kernel_sum.setArg(18, buf_phase);
    kernel_sum.setArg(19, buf_stress_array);
    kernel_sum.setArg(20, buf_wall_velocity);

    kernel_distances.setArg(0, buf_x0);
    kernel_distances.setArg(1, buf_y0);
    kernel_distances.setArg(2, buf_z0);
    kernel_distances.setArg(3, buf_particle_0);
    kernel_distances.setArg(4, buf_particle_1);
    kernel_distances.setArg(5, buf_length);
    kernel_distances.setArg(6, buf_valid);
    kernel_distances.setArg(7, buf_dimensions);

    kernel_validate.setArg(0, buf_valid);
    kernel_validate.setArg(1, buf_x0);
    kernel_validate.setArg(2, buf_y0);
    kernel_validate.setArg(3, buf_z0);
    kernel_validate.setArg(4, buf_x1);
    kernel_validate.setArg(5, buf_y1);
    kernel_validate.setArg(6, buf_z1);
    kernel_validate.setArg(7, buf_original_delta_t);
    kernel_validate.setArg(8, buf_delta_t);
    kernel_validate.setArg(9, buf_t);
    kernel_validate.setArg(10, buf_magnetic_field);
    kernel_validate.setArg(11, buf_mason);
    kernel_validate.setArg(12, buf_ar);
    kernel_validate.setArg(13, buf_dimensions);
    kernel_validate.setArg(14, buf_mode);
    kernel_validate.setArg(15, buf_phase);
    kernel_validate.setArg(16, buf_stress_array);
    kernel_validate.setArg(17, buf_stress);
    kernel_validate.setArg(18, buf_particles);
    kernel_validate.setArg(19, buf_length);
    kernel_validate.setArg(20, buf_field_direction);
}

void SimulationContextOCL::setPositions(
    const double* x,
    const double* y,
    const double* z)
{
    m_queue.enqueueWriteBuffer(buf_x0, CL_TRUE, 0, m_particles * sizeof(double), x);
    m_queue.enqueueWriteBuffer(buf_y0, CL_TRUE, 0, m_particles * sizeof(double), y);
    m_queue.enqueueWriteBuffer(buf_z0, CL_TRUE, 0, m_particles * sizeof(double), z);
}

void SimulationContextOCL::readPositions(
    double* x,
    double* y,
    double* z)
{
    m_queue.enqueueReadBuffer(buf_x1, CL_TRUE, 0, m_particles * sizeof(double), x);
    m_queue.enqueueReadBuffer(buf_y1, CL_TRUE, 0, m_particles * sizeof(double), y);
    m_queue.enqueueReadBuffer(buf_z1, CL_TRUE, 0, m_particles * sizeof(double), z);
}

void SimulationContextOCL::setPhase(int phase){
    m_queue.enqueueWriteBuffer(buf_phase, CL_TRUE, 0, sizeof(int), &phase);
}
void SimulationContextOCL::setMode(int mode) {
    m_queue.enqueueWriteBuffer(buf_mode, CL_TRUE, 0, sizeof(int), &mode);
}
void SimulationContextOCL::setMagneticField(double magnetic_field[3]) {
    m_queue.enqueueWriteBuffer(buf_magnetic_field, CL_TRUE, 0, 3*sizeof(double), &magnetic_field);
}
void SimulationContextOCL::setMason(double mason){
    m_queue.enqueueWriteBuffer(buf_mason, CL_TRUE, 0, sizeof(double), &mason);
}
void SimulationContextOCL::setWallVelocity(double wall_velocity){
    m_queue.enqueueWriteBuffer(buf_wall_velocity, CL_TRUE, 0, sizeof(double), &wall_velocity);
}
void SimulationContextOCL::readValid(int* valid){
    m_queue.enqueueReadBuffer(buf_valid, CL_TRUE, 0, sizeof(int), valid);
}
void SimulationContextOCL::readTime(double* time){
    m_queue.enqueueReadBuffer(buf_t, CL_TRUE, 0, sizeof(double), time);
}
void SimulationContextOCL::readDeltaT(double* delta_t){
    m_queue.enqueueReadBuffer(buf_delta_t, CL_TRUE, 0, sizeof(double), delta_t);
}
void SimulationContextOCL::readStress(double* stress){
    m_queue.enqueueReadBuffer(buf_stress, CL_TRUE, 0, sizeof(double), stress);
}

void SimulationContextOCL::enqueueStep() {
    m_queue.enqueueNDRangeKernel(
        kernel_validate,
        cl::NullRange,
        cl::NDRange(m_particles),
        cl::NullRange
    );

    m_queue.enqueueNDRangeKernel(
        kernel_forces,
        cl::NullRange,
        cl::NDRange(m_matrix_size),
        cl::NullRange
    );

    m_queue.enqueueNDRangeKernel(
        kernel_sum,
        cl::NullRange,
        cl::NDRange(m_particles),
        cl::NullRange
    );

    m_queue.enqueueNDRangeKernel(
        kernel_distances,
        cl::NullRange,
        cl::NDRange(m_matrix_size),
        cl::NullRange
    );
}