#pragma once
#include <CL/cl.hpp>
#include "SubmitJob.h"
#include "SubmitQueue.h"


class SimulationContextOCL {
public:
    SimulationContextOCL(int particles, int dimensions, int length, double field_direction, double ma, double ar, double original_delta_t);

    void setPositions(const double* x, const double* y, const double* z);
    void readPositions(double* x, double* y, double* z);
    void setPhase(int phase);
    void setMode(int mode);
    void setMagneticField(double magnetic_field[3]);
    void setMason(double mason);
    void setWallVelocity(double wall_velocity);
    int readValid();
    double readTime();
    void readDeltaT(double* delta_t);
    void readStress(double* stress);

    StepSync enqueueStep();

private:
    void createBuffers();
    void loadKernels();

    int m_particles;
    int m_dimensions;
    int m_length;
    int m_matrix_size;
    int m_valid;
    double m_mason;
    double m_ar;
    double m_time;
    double m_delta_t;
    double m_field_direction;

    cl::Context& m_context;
    cl::Device& m_device;
    cl::CommandQueue& m_queue;

    cl::Buffer buf_x0, buf_y0, buf_z0;
    cl::Buffer buf_x1, buf_y1, buf_z1;
    cl::Buffer buf_fx, buf_fy, buf_fz;

    cl::Buffer buf_magnetic_field, buf_field_direction, buf_wall_velocity;
    cl::Buffer buf_mode, buf_length, buf_phase;
    cl::Buffer buf_dimensions, buf_particles, buf_matrix_size;
    cl::Buffer buf_t, buf_delta_t, buf_original_delta_t;
    cl::Buffer buf_valid, buf_mason, buf_ar;
    cl::Buffer buf_initial_idx, buf_last_idx, buf_particle_0;
    cl::Buffer buf_stress, buf_stress_array, buf_particle_1;

    cl::Program prog_forces;
    cl::Program prog_sum;
    cl::Program prog_validate;
    cl::Program prog_distances;
    cl::Kernel  kernel_forces;
    cl::Kernel  kernel_sum;
    cl::Kernel  kernel_validate;
    cl::Kernel  kernel_distances;
};
