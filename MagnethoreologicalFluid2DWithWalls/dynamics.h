using namespace std::chrono;
void Simulation(int particles, int length, double mason, double amplitude_relationship, double original_delta_t, int repetition, double max_time) {
	auto start = high_resolution_clock::now();
	std::cout << "Starting simulation for AR = " + std::to_string(amplitude_relationship) + " and Mason = " + std::to_string(mason) + "\n";

	double magnetic_field[2] = { 0.0, 1.0 }; 
	double delta_t = original_delta_t;
	double pi = 3.14159265359;
	double step = 2 * pi / (mason * 360);
	int laps = ceil(max_time * mason / (2 * pi));
	double time = 0.0;
	int lap = 0;
	int current_lap = 0;
	int counter = 0;
	int window = 5;
	int matrix_size = particles * (particles - 1) / 2;
	int valid = 1;
	int* initial_indices_sum = new int[particles];
	int* last_indices_sum = new int[particles];
	int* particle_0 = new int[matrix_size];
	int* particle_1 = new int[matrix_size];
	bool end_simulation = false;

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
			particle_0[index] = i;
			particle_1[index] = j;
			index++;
		}
	}

	Analysis analysis(mason, amplitude_relationship, particles, length, window);
	Box box(particles, length);
	box.WritePositions(counter, mason, amplitude_relationship, repetition, "");

	std::vector<double> get_x = box.ReturnX();
	double* x_0 = get_x.data();
	std::vector<double> get_y = box.ReturnY();
	double* y_0 = get_y.data();

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

	cl::Buffer buffer_x_0 = cl::Buffer(context, CL_MEM_READ_WRITE, particles * sizeof(double));
	cl::Buffer buffer_y_0 = cl::Buffer(context, CL_MEM_READ_WRITE, particles * sizeof(double));
	cl::Buffer buffer_magnetic_field = cl::Buffer(context, CL_MEM_READ_WRITE, 2 * sizeof(double));
	cl::Buffer buffer_length = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_particles = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_delta_t = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double));
	cl::Buffer buffer_original_delta_t = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(double));
	cl::Buffer buffer_time = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double));
	cl::Buffer buffer_forces_x = cl::Buffer(context, CL_MEM_READ_WRITE, matrix_size * sizeof(double));
	cl::Buffer buffer_forces_y = cl::Buffer(context, CL_MEM_READ_WRITE, matrix_size * sizeof(double));
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
	cl::Buffer buffer_r_array = cl::Buffer(context, CL_MEM_READ_WRITE, matrix_size * sizeof(double));

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

	sum_kernel.setArg(0, buffer_x_0);
	sum_kernel.setArg(1, buffer_y_0);
	sum_kernel.setArg(2, buffer_length);
	sum_kernel.setArg(3, buffer_particles);
	sum_kernel.setArg(4, buffer_delta_t);
	sum_kernel.setArg(5, buffer_forces_x);
	sum_kernel.setArg(6, buffer_forces_y);
	sum_kernel.setArg(7, buffer_initial_indices_sum);
	sum_kernel.setArg(8, buffer_last_indices_sum);
	sum_kernel.setArg(9, buffer_matrix_size);
	sum_kernel.setArg(10, buffer_valid);
	sum_kernel.setArg(11, buffer_x_1);
	sum_kernel.setArg(12, buffer_y_1);

	forces_kernel.setArg(0, buffer_x_0);
	forces_kernel.setArg(1, buffer_y_0);
	forces_kernel.setArg(2, buffer_magnetic_field);
	forces_kernel.setArg(3, buffer_length);
	forces_kernel.setArg(4, buffer_particles);
	forces_kernel.setArg(5, buffer_matrix_size);
	forces_kernel.setArg(6, buffer_particle_0);
	forces_kernel.setArg(7, buffer_particle_1);
	forces_kernel.setArg(8, buffer_forces_x);
	forces_kernel.setArg(9, buffer_forces_y);

	distances_kernel.setArg(0, buffer_x_0);
	distances_kernel.setArg(1, buffer_y_0);
	distances_kernel.setArg(2, buffer_particle_0);
	distances_kernel.setArg(3, buffer_particle_1);
	distances_kernel.setArg(4, buffer_length);
	distances_kernel.setArg(5, buffer_valid);

	validation_kernel.setArg(0, buffer_valid);
	validation_kernel.setArg(1, buffer_x_0);
	validation_kernel.setArg(2, buffer_y_0);
	validation_kernel.setArg(3, buffer_x_1);
	validation_kernel.setArg(4, buffer_y_1);
	validation_kernel.setArg(5, buffer_original_delta_t);
	validation_kernel.setArg(6, buffer_delta_t);
	validation_kernel.setArg(7, buffer_time);
	validation_kernel.setArg(8, buffer_magnetic_field);
	validation_kernel.setArg(9, buffer_mason);
	validation_kernel.setArg(10, buffer_amplitude_relationship);

	cl::NDRange global_long(matrix_size);
	cl::NDRange global_short(particles);

	queue.enqueueWriteBuffer(buffer_particles, CL_TRUE, 0, sizeof(int), &particles);
	queue.enqueueWriteBuffer(buffer_length, CL_TRUE, 0, sizeof(int), &length);
	queue.enqueueWriteBuffer(buffer_matrix_size, CL_TRUE, 0, sizeof(int), &matrix_size);
	queue.enqueueWriteBuffer(buffer_delta_t, CL_TRUE, 0, sizeof(double), &delta_t);
	queue.enqueueWriteBuffer(buffer_original_delta_t, CL_TRUE, 0, sizeof(double), &delta_t);
	queue.enqueueWriteBuffer(buffer_mason, CL_TRUE, 0, sizeof(double), &mason);
	queue.enqueueWriteBuffer(buffer_amplitude_relationship, CL_TRUE, 0, sizeof(double), &amplitude_relationship);
	queue.enqueueWriteBuffer(buffer_time, CL_TRUE, 0, sizeof(double), &time);
	queue.enqueueWriteBuffer(buffer_magnetic_field, CL_TRUE, 0, 2 * sizeof(double), &magnetic_field);
	queue.enqueueWriteBuffer(buffer_particle_0, CL_TRUE, 0, matrix_size * sizeof(int), particle_0);
	queue.enqueueWriteBuffer(buffer_particle_1, CL_TRUE, 0, matrix_size * sizeof(int), particle_1);
	queue.enqueueWriteBuffer(buffer_initial_indices_sum, CL_TRUE, 0, particles * sizeof(int), initial_indices_sum);
	queue.enqueueWriteBuffer(buffer_last_indices_sum, CL_TRUE, 0, particles * sizeof(int), last_indices_sum);
	queue.enqueueWriteBuffer(buffer_x_0, CL_TRUE, 0, particles * sizeof(double), x_0);
	queue.enqueueWriteBuffer(buffer_y_0, CL_TRUE, 0, particles * sizeof(double), y_0);

	while (!end_simulation) {
		queue.enqueueNDRangeKernel(forces_kernel, cl::NullRange, global_long, cl::NullRange);
		queue.enqueueNDRangeKernel(sum_kernel, cl::NullRange, global_short, cl::NullRange);
		queue.enqueueNDRangeKernel(distances_kernel, cl::NullRange, global_long, cl::NullRange);
		queue.enqueueNDRangeKernel(validation_kernel, cl::NullRange, global_short, cl::NullRange);
	
		queue.enqueueReadBuffer(buffer_time, CL_TRUE, 0, sizeof(double), &time);

		current_lap = floor(time * mason / (2 * pi));
		if (current_lap > lap) {
			counter++;
			lap = current_lap;
			queue.enqueueReadBuffer(buffer_x_0, CL_TRUE, 0, particles * sizeof(double), x_0);
			queue.enqueueReadBuffer(buffer_y_0, CL_TRUE, 0, particles * sizeof(double), y_0);
			if (repetition == 0) {
				box.SetX(x_0);
				box.SetY(y_0);
				box.WritePositions(counter, mason, amplitude_relationship, repetition, "");
			}
			end_simulation = analysis.PreAnalysis(x_0, y_0, time) || (current_lap > laps);
		}
	}
	analysis.WriteAnalysis(repetition, "");

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	std::cout << "Finishing simulation for thread " << repetition << ". Took " << duration.count() << " seconds." << "\n";
}