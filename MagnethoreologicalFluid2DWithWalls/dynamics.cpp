#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 200

#include "dynamics.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <chrono>
#include <random>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <CL/opencl.hpp>
#include <thread>
#include <cmath>
#include "box.h"
#include "analysis.h"

using namespace std::chrono;
using namespace cl;
void Simulation(double field_direction, int phases, int particles, int dimensions, int length, double mason, 
	double amplitude_relationship, double original_delta_t, int repetition, double max_times[3],
	bool keep_positions, bool load_positions, double creep_time) {
	auto start = high_resolution_clock::now();
	if (load_positions) {
		std::cout << "Starting simulation for AR = " + std::to_string(amplitude_relationship) + ", Mason = " + std::to_string(mason) + ", field direction = " + std::to_string(field_direction) + ", creep time " + std::to_string(creep_time) + " and repetition " + std::to_string(repetition) + "\n";
	}
	else {
		std::cout << "Starting simulation for AR = " + std::to_string(amplitude_relationship) + ", Mason = " + std::to_string(mason) + ", field direction = " + std::to_string(field_direction) + " and repetition " + std::to_string(repetition) + "\n";
	}

	double magnetic_field[3];
	if (dimensions == 3) {
		magnetic_field[0] = 0.0;
		magnetic_field[1] = 0.0;
		magnetic_field[2] = 1.0;
	}
	else {
		magnetic_field[0] = 0.0;
		magnetic_field[1] = 1.0;
		magnetic_field[2] = 0.0;
	}

	double frecuency = mason + (mason < 0.00000001) * 1.0;
	double delta_t = original_delta_t;
	double pi = 3.14159265359;
	//step divides a whole cycle depending on the frequency (mason number).
	double step = 2 * pi / (mason * 360);
	//strech tracks if the simulation overcame a step. This is used during the last phase of the simulation to keep track of the changes while applying stress.
	double stretch = 0;
	//Max time for a phase.
	double max_time;
	//Current time of the simulation.
	double time = max_times[phases - 2] * load_positions;
	//Variable to keep track of the time during stress phase.
	double t = 0;
	//Number laps during a phase.
	int laps;
	//Variable to check if current_lap changed. Logic behind this: if current_lap > lap, then, we are in a new lap.
	int lap = 0;
	//Current lap of the phase.
	int current_lap = 0;
	//Counter to record stresses.
	int counter = 0;
	int window = 5;
	int matrix_size = particles * (particles - 1) / 2;
	int valid = 1;
	int mode = 0;
	int* initial_indices_sum = new int[particles];
	int* last_indices_sum = new int[particles];
	int* particle_0 = new int[matrix_size];
	int* particle_1 = new int[matrix_size];
	bool end_simulation;
	double stress = 0;
	int file_to_load = phases == 3 ? ceil(max_times[0] * frecuency / (2 * pi)) + ceil(max_times[1] * frecuency / (2 * pi)) : ceil(max_times[0] * frecuency / (2 * pi));
	double wall_velocity = length;
	//These variables are meant for creep experiment. First one is to set if we are in relaxation time, 0, or not, 1.
	int relaxation = 1;
	//To keep track of how much time has passed during the steps while creeping.
	double creep_chrono = 0.0;

	for (int i = 0; i < particles; i++) {
		initial_indices_sum[i] = 0;
		last_indices_sum[i] = 0;
	}
	for (int i = 1; i < particles; i++) {
		for (int j = 1; j <= i; j++) {
			initial_indices_sum[i] += particles - j;
		}
	}
	last_indices_sum[0] = particles - 1;
	for (int i = 1; i < particles; i++) {
		for (int j = 1; j <= i + 1; j++) {
			last_indices_sum[i] += particles - j;
		}
	}

	int index = 0;
	for (int i = 0; i < particles - 1; i++) {
		for (int j = i + 1; j < particles; j++) {
			particle_1[index] = i;
			particle_0[index] = j;
			index++;
		}
	}

	std::string tag = load_positions ? "field_direction-" + std::to_string(field_direction) + "-creep_time-" + std::to_string(creep_time) : "field_direction-" + std::to_string(field_direction);
	Analysis* analysis = new Analysis(mason, amplitude_relationship, particles, length, window, dimensions, field_direction);
	Box* box = new Box(particles, length, dimensions);
	if (load_positions) box->ReadCsv("positions/positions-" + std::to_string(mason) + "-" + std::to_string(amplitude_relationship) + "-" + std::to_string(repetition) + "-" + std::to_string(file_to_load) + "-field_direction-" + std::to_string(field_direction) + ".csv");

	std::vector<double> get_x = box->ReturnX();
	std::vector<double> get_y = box->ReturnY();
	std::vector<double> get_z = box->ReturnZ();

	double r_min = 1 - log(100.0) / 10;

	std::vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);
	if (all_platforms.size() == 0) {
		std::cout << " No platform found. Check OpenCL installation.\n";
		exit(1);
	}
	cl::Platform platform = all_platforms[0];
	std::vector<cl::Device> all_devices;
	platform.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);
	if (all_devices.size() == 0) {
		std::cout << "No devices found. Check OpenCL installation.\n";
		exit(1);
	}
	cl::Device device = all_devices[0];
	cl::Context context({ device });
	cl::CommandQueue queue = cl::CommandQueue(context, device);

	// SVM allocations
	cl::SVMAllocator<double, cl::SVMTraitFine<>> svmAlloc(context);
	double* x_0 = svmAlloc.allocate(particles);
	double* y_0 = svmAlloc.allocate(particles);
	double* z_0 = svmAlloc.allocate(particles);
	for (int i = 0; i < particles; i++) {
		x_0[i] = get_x.data()[i];
		y_0[i] = get_y.data()[i];
		z_0[i] = get_z.data()[i];
	}


	cl::Buffer buffer_magnetic_field = cl::Buffer(context, CL_MEM_READ_WRITE, 3 * sizeof(double));
	cl::Buffer buffer_mode = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_phase = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_dimensions = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_length = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_field_direction = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(double));
	cl::Buffer buffer_particles = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_delta_t = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double));
	cl::Buffer buffer_original_delta_t = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(double));
	cl::Buffer buffer_time = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double));
	cl::Buffer buffer_forces_x = cl::Buffer(context, CL_MEM_READ_WRITE, matrix_size * sizeof(double));
	cl::Buffer buffer_forces_y = cl::Buffer(context, CL_MEM_READ_WRITE, matrix_size * sizeof(double));
	cl::Buffer buffer_forces_z = cl::Buffer(context, CL_MEM_READ_WRITE, matrix_size * sizeof(double));
	cl::Buffer buffer_initial_indices_sum = cl::Buffer(context, CL_MEM_READ_ONLY, matrix_size * sizeof(int));
	cl::Buffer buffer_last_indices_sum = cl::Buffer(context, CL_MEM_READ_ONLY, matrix_size * sizeof(int));
	cl::Buffer buffer_particle_0 = cl::Buffer(context, CL_MEM_READ_ONLY, matrix_size * sizeof(int));
	cl::Buffer buffer_particle_1 = cl::Buffer(context, CL_MEM_READ_ONLY, matrix_size * sizeof(int));
	cl::Buffer buffer_matrix_size = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_valid = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(int));
	cl::Buffer buffer_mason = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(double));
	cl::Buffer buffer_amplitude_relationship = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(double));
	cl::Buffer buffer_x_1 = cl::Buffer(context, CL_MEM_READ_WRITE, particles * sizeof(double));
	cl::Buffer buffer_y_1 = cl::Buffer(context, CL_MEM_READ_WRITE, particles * sizeof(double));
	cl::Buffer buffer_z_1 = cl::Buffer(context, CL_MEM_READ_WRITE, particles * sizeof(double));
	cl::Buffer buffer_stress_array = cl::Buffer(context, CL_MEM_READ_WRITE, particles * sizeof(double));
	cl::Buffer buffer_stress = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double));
	cl::Buffer buffer_r_array = cl::Buffer(context, CL_MEM_READ_WRITE, matrix_size * sizeof(double));
	cl::Buffer buffer_wall_velocity = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(double));


	std::ifstream sum_file("sum.cl");
	std::string sum_code(std::istreambuf_iterator<char>(sum_file), (std::istreambuf_iterator<char>()));
	cl::Program::Sources sum_source;
	sum_source.push_back({ sum_code.c_str(), sum_code.length() });

	std::ifstream forces_file("forces.cl");
	std::string forces_code(std::istreambuf_iterator<char>(forces_file), (std::istreambuf_iterator<char>()));
	cl::Program::Sources forces_source;
	forces_source.push_back({ forces_code.c_str(), forces_code.length() });

	std::ifstream distances_file("distances.cl");
	std::string distances_code(std::istreambuf_iterator<char>(distances_file), (std::istreambuf_iterator<char>()));
	cl::Program::Sources distances_source;
	distances_source.push_back({ distances_code.c_str(), distances_code.length() });

	std::ifstream validation_file("validation.cl");
	std::string validation_code(std::istreambuf_iterator<char>(validation_file), (std::istreambuf_iterator<char>()));
	cl::Program::Sources validation_source;
	validation_source.push_back({ validation_code.c_str(), validation_code.length() });

	cl::Program sum_program = cl::Program(context, sum_source);
	if (sum_program.build({ device }) != CL_SUCCESS) {
		std::string problem = sum_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
		std::cout << "Error building: " << problem << "\n";
		exit(1);
	}

	cl::Program forces_program = cl::Program(context, forces_source);
	if (forces_program.build({ device }) != CL_SUCCESS) {
		std::string problem = forces_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
		std::cout << "Error building: " << problem << "\n";
		exit(1);
	}

	cl::Program distances_program = cl::Program(context, distances_source);
	if (distances_program.build({ device }) != CL_SUCCESS) {
		std::string problem = distances_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
		std::cout << "Error building: " << problem << "\n";
		exit(1);
	}

	cl::Program validation_program = cl::Program(context, validation_source);
	if (validation_program.build({ device }) != CL_SUCCESS) {
		std::string problem = validation_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
		std::cout << "Error building: " << problem << "\n";
		exit(1);
	}

	cl::Kernel sum_kernel(sum_program, "sum");
	cl::Kernel forces_kernel(forces_program, "forces");
	cl::Kernel distances_kernel(distances_program, "distances");
	cl::Kernel validation_kernel(validation_program, "validation");

	sum_kernel.setArg(0, x_0);
	sum_kernel.setArg(1, y_0);
	sum_kernel.setArg(2, z_0);
	sum_kernel.setArg(3, buffer_dimensions);
	sum_kernel.setArg(4, buffer_length);
	sum_kernel.setArg(5, buffer_particles);
	sum_kernel.setArg(6, buffer_delta_t);
	sum_kernel.setArg(7, buffer_forces_x);
	sum_kernel.setArg(8, buffer_forces_y);
	sum_kernel.setArg(9, buffer_forces_z);
	sum_kernel.setArg(10, buffer_initial_indices_sum);
	sum_kernel.setArg(11, buffer_last_indices_sum);
	sum_kernel.setArg(12, buffer_matrix_size);
	sum_kernel.setArg(13, buffer_valid);
	sum_kernel.setArg(14, buffer_x_1);
	sum_kernel.setArg(15, buffer_y_1);
	sum_kernel.setArg(16, buffer_z_1);
	sum_kernel.setArg(17, buffer_mode);
	sum_kernel.setArg(18, buffer_phase);
	sum_kernel.setArg(19, buffer_stress_array);
	sum_kernel.setArg(20, buffer_wall_velocity);

	forces_kernel.setArg(0, x_0);
	forces_kernel.setArg(1, y_0);
	forces_kernel.setArg(2, z_0);
	forces_kernel.setArg(3, buffer_dimensions);
	forces_kernel.setArg(4, buffer_magnetic_field);
	forces_kernel.setArg(5, buffer_length);
	forces_kernel.setArg(6, buffer_particles);
	forces_kernel.setArg(7, buffer_matrix_size);
	forces_kernel.setArg(8, buffer_particle_0);
	forces_kernel.setArg(9, buffer_particle_1);
	forces_kernel.setArg(10, buffer_forces_x);
	forces_kernel.setArg(11, buffer_forces_y);
	forces_kernel.setArg(12, buffer_forces_z);

	distances_kernel.setArg(0, x_0);
	distances_kernel.setArg(1, y_0);
	distances_kernel.setArg(2, z_0);
	distances_kernel.setArg(3, buffer_particle_0);
	distances_kernel.setArg(4, buffer_particle_1);
	distances_kernel.setArg(5, buffer_length);
	distances_kernel.setArg(6, buffer_valid);
	distances_kernel.setArg(7, buffer_dimensions);

	validation_kernel.setArg(0, buffer_valid);
	validation_kernel.setArg(1, x_0);
	validation_kernel.setArg(2, y_0);
	validation_kernel.setArg(3, z_0);
	validation_kernel.setArg(4, buffer_x_1);
	validation_kernel.setArg(5, buffer_y_1);
	validation_kernel.setArg(6, buffer_z_1);
	validation_kernel.setArg(7, buffer_original_delta_t);
	validation_kernel.setArg(8, buffer_delta_t);
	validation_kernel.setArg(9, buffer_time);
	validation_kernel.setArg(10, buffer_magnetic_field);
	validation_kernel.setArg(11, buffer_mason);
	validation_kernel.setArg(12, buffer_amplitude_relationship);
	validation_kernel.setArg(13, buffer_dimensions);
	validation_kernel.setArg(14, buffer_mode);
	validation_kernel.setArg(15, buffer_phase);
	validation_kernel.setArg(16, buffer_stress_array);
	validation_kernel.setArg(17, buffer_stress);
	validation_kernel.setArg(18, buffer_particles);
	validation_kernel.setArg(19, buffer_length);
	validation_kernel.setArg(20, buffer_field_direction);

	cl::NDRange global_long(matrix_size);
	cl::NDRange global_short(particles);

	queue.enqueueWriteBuffer(buffer_particles, CL_TRUE, 0, sizeof(int), &particles);
	queue.enqueueWriteBuffer(buffer_dimensions, CL_TRUE, 0, sizeof(int), &dimensions);
	queue.enqueueWriteBuffer(buffer_length, CL_TRUE, 0, sizeof(int), &length);
	queue.enqueueWriteBuffer(buffer_field_direction, CL_TRUE, 0, sizeof(double), &field_direction);
	queue.enqueueWriteBuffer(buffer_matrix_size, CL_TRUE, 0, sizeof(int), &matrix_size);
	queue.enqueueWriteBuffer(buffer_delta_t, CL_TRUE, 0, sizeof(double), &delta_t);
	queue.enqueueWriteBuffer(buffer_original_delta_t, CL_TRUE, 0, sizeof(double), &delta_t);
	queue.enqueueWriteBuffer(buffer_original_delta_t, CL_TRUE, 0, sizeof(double), &original_delta_t);
	queue.enqueueWriteBuffer(buffer_mason, CL_TRUE, 0, sizeof(double), &mason);
	queue.enqueueWriteBuffer(buffer_amplitude_relationship, CL_TRUE, 0, sizeof(double), &amplitude_relationship);
	queue.enqueueWriteBuffer(buffer_time, CL_TRUE, 0, sizeof(double), &time);
	queue.enqueueWriteBuffer(buffer_magnetic_field, CL_TRUE, 0, 3 * sizeof(double), &magnetic_field);
	queue.enqueueWriteBuffer(buffer_particle_0, CL_TRUE, 0, matrix_size * sizeof(int), particle_0);
	queue.enqueueWriteBuffer(buffer_particle_1, CL_TRUE, 0, matrix_size * sizeof(int), particle_1);
	queue.enqueueWriteBuffer(buffer_initial_indices_sum, CL_TRUE, 0, particles * sizeof(int), initial_indices_sum);
	queue.enqueueWriteBuffer(buffer_last_indices_sum, CL_TRUE, 0, particles * sizeof(int), last_indices_sum);

	for (int phase = 0 + (phases - 1)*load_positions; phase < phases; phase++) {
		max_time = max_times[phase];
		laps = ceil(max_time * frecuency / (2 * pi));
		queue.enqueueWriteBuffer(buffer_phase, CL_TRUE, 0, sizeof(int), &phase);
		mode = phase == phases - 1;
		queue.enqueueWriteBuffer(buffer_mode, CL_TRUE, 0, sizeof(int), &mode);
		end_simulation = delta_t < original_delta_t / 16.0;

		if (phases > 2 && phase == 1) {
			queue.enqueueWriteBuffer(buffer_magnetic_field, CL_TRUE, 0, 3 * sizeof(double), &magnetic_field);
			double perturbation = 0.0;
			queue.enqueueWriteBuffer(buffer_mason, CL_TRUE, 0, sizeof(double), &perturbation);
			queue.enqueueWriteBuffer(buffer_wall_velocity, CL_TRUE, 0, sizeof(double), &wall_velocity);
		}

		while (!end_simulation) {
			queue.enqueueNDRangeKernel(forces_kernel, cl::NullRange, global_long, cl::NullRange);
			queue.enqueueNDRangeKernel(sum_kernel, cl::NullRange, global_short, cl::NullRange);
			queue.enqueueNDRangeKernel(distances_kernel, cl::NullRange, global_long, cl::NullRange);
			queue.enqueueReadBuffer(buffer_valid, CL_TRUE, 0, sizeof(int), &valid);
			queue.enqueueNDRangeKernel(validation_kernel, cl::NullRange, global_short, cl::NullRange);
			queue.enqueueReadBuffer(buffer_time, CL_TRUE, 0, sizeof(double), &time);

			current_lap = floor(time * frecuency / (2 * pi));
			if (current_lap > lap) {
				counter++;
				lap = current_lap;
				if (keep_positions) {
					box->SetX(x_0);
					box->SetY(y_0);
					box->SetZ(z_0);
					box->WritePositions(counter, mason, amplitude_relationship, repetition, tag);
				}
				analysis->PreAnalysis(x_0, y_0, z_0, time);
				analysis->Connectivity(x_0, y_0, z_0);
				end_simulation = time > max_time;
			}
			else if (current_lap == laps - 1) {
				queue.enqueueReadBuffer(buffer_delta_t, CL_TRUE, 0, sizeof(double), &delta_t);
				stretch += delta_t;

				if (stretch > step) {
					counter++;
					if (keep_positions) {
						box->SetX(x_0);
						box->SetY(y_0);
						box->SetZ(z_0);
						box->WritePositions(counter, mason, amplitude_relationship, repetition, tag);
					}
					analysis->PreAnalysis(x_0, y_0, z_0, time);
					analysis->Connectivity(x_0, y_0, z_0);
					end_simulation = time > max_time;
					stretch = 0;
				}
			}
			if (phase == phases - 1 && valid == 1) {
				t += delta_t;
				creep_chrono += delta_t;
				//If we are running creep experiment and passed time is equal or greater than creep phase, change relaxation status.
				relaxation = !(creep_chrono >= creep_time && load_positions);

				queue.enqueueReadBuffer(buffer_stress, CL_TRUE, 0, sizeof(double), &stress);
				analysis->RecordStress(t, stress);
				end_simulation = time > max_time;
				wall_velocity = length*relaxation;
				queue.enqueueWriteBuffer(buffer_wall_velocity, CL_TRUE, 0, sizeof(double), &wall_velocity);
			}

			if (delta_t < original_delta_t / 16.0) {
				end_simulation = true;

				std::ofstream file{ "failed_simulations.txt" };
				file << "Failed simulation for Ma " + std::to_string(mason) + ", AR " + std::to_string(amplitude_relationship) + ", field direction = " + std::to_string(field_direction) + ", creep time " + std::to_string(creep_time) + " and repetition " + std::to_string(repetition) + "\n";
			}
		}
	}

	if (keep_positions) {
		box->SetX(x_0);
		box->SetY(y_0);
		box->SetZ(z_0);
		box->WritePositions(counter, mason, amplitude_relationship, repetition, tag);
	}
	analysis->WriteAnalysis(repetition, tag);
	analysis->WriteStress(repetition, tag);
	analysis->WriteConnectivity(repetition, tag);

	delete box;
	delete analysis;
	delete[] initial_indices_sum;
	delete[] last_indices_sum;
	delete[] particle_0;
	delete[] particle_1;

	svmAlloc.deallocate(x_0, particles);
	svmAlloc.deallocate(y_0, particles);
	svmAlloc.deallocate(z_0, particles);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	std::cout << "Finishing simulation for AR = " + std::to_string(amplitude_relationship) + ", Mason = " + std::to_string(mason) + ", field direction " + std::to_string(field_direction) + " and repetition " + std::to_string(repetition) + ". Took " + std::to_string(duration.count()) + " seconds.\n";
}