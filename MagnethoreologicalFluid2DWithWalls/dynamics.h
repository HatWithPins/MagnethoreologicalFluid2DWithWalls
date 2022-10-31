void TestSum() {
	double x_0[4] = { 0, 0, 0, 0 };
	double x_1[4];
	double y_0[4] = { 0, 0, 0, 0 };
	double y_1[4];
	int length = 10;
	int particles = 4;
	double delta_t = 1.0;
	double forces_x[6] = { 1, 2, 3, 4, 5, 6 };
	double forces_y[6] = { 1, 2, 3, 4, 5, 6 };
	int matrix_size = 6;
	int valid = 1;
	int initial_indices_sum[4];
	int last_indices_sum[4];

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

	int expected_initial[4] = { 0, 3, 5, 6 };
	int expected_last[4] = { 3, 5, 6, 6 };
	for (int i = 0; i < particles; i++) {
		std::cout << "Expected initial index: " << expected_initial[i] << " Actual index: " << initial_indices_sum[i] << "\n";
		std::cout << "Expected last index: " << expected_last[i] << " Actual index: " << last_indices_sum[i] << "\n";
	}

	double expected_sum[4] = { 6, 8, 0, -14 };

	std::vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);
	if (all_platforms.size() == 0) {
		std::cout << " No se encontró ninguna plataforma. Comprueba la instalación de OpenCL.\n";
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
	cl::Buffer buffer_length = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_particles = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_delta_t = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(double));
	cl::Buffer buffer_forces_x = cl::Buffer(context, CL_MEM_READ_WRITE, matrix_size * sizeof(double));
	cl::Buffer buffer_forces_y = cl::Buffer(context, CL_MEM_READ_WRITE, matrix_size * sizeof(double));
	cl::Buffer buffer_initial_indices_sum = cl::Buffer(context, CL_MEM_READ_ONLY, matrix_size * sizeof(int));
	cl::Buffer buffer_last_indices_sum = cl::Buffer(context, CL_MEM_READ_ONLY, matrix_size * sizeof(int));
	cl::Buffer buffer_matrix_size = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int));
	cl::Buffer buffer_valid = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(int));
	cl::Buffer buffer_x_1 = cl::Buffer(context, CL_MEM_READ_WRITE, particles * sizeof(double));
	cl::Buffer buffer_y_1 = cl::Buffer(context, CL_MEM_READ_WRITE, particles * sizeof(double));

	std::ifstream sum_file("sum.cl");
	std::string sum_code(std::istreambuf_iterator<char>(sum_file), (std::istreambuf_iterator<char>()));
	cl::Program::Sources sum_source;
	sum_source.push_back({ sum_code.c_str(), sum_code.length() });

	cl::Program sum_program = cl::Program(context, sum_source);
	if (sum_program.build({ device }) != CL_SUCCESS) {
		std::string problem = sum_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
		std::cout << "Error building: " << problem << "\n";
		exit(1);
	}

	cl::Kernel sum_kernel(sum_program, "sum");

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

	queue.enqueueWriteBuffer(buffer_x_0, CL_TRUE, 0, particles * sizeof(double), &x_0);
	queue.enqueueWriteBuffer(buffer_y_0, CL_TRUE, 0, particles * sizeof(double), &y_0);
	queue.enqueueWriteBuffer(buffer_forces_x, CL_TRUE, 0, matrix_size * sizeof(double), &forces_x);
	queue.enqueueWriteBuffer(buffer_forces_y, CL_TRUE, 0, matrix_size * sizeof(double), &forces_y);
	queue.enqueueWriteBuffer(buffer_initial_indices_sum, CL_TRUE, 0, particles * sizeof(int), &initial_indices_sum);
	queue.enqueueWriteBuffer(buffer_last_indices_sum, CL_TRUE, 0, particles * sizeof(int), &last_indices_sum);
	queue.enqueueWriteBuffer(buffer_length, CL_TRUE, 0, sizeof(int), &length);
	queue.enqueueWriteBuffer(buffer_particles, CL_TRUE, 0, sizeof(int), &particles);
	queue.enqueueWriteBuffer(buffer_matrix_size, CL_TRUE, 0, sizeof(int), &matrix_size);
	queue.enqueueWriteBuffer(buffer_valid, CL_TRUE, 0, sizeof(int), &valid);
	queue.enqueueWriteBuffer(buffer_delta_t, CL_TRUE, 0, sizeof(double), &delta_t);

	cl::NDRange global_long(matrix_size);
	cl::NDRange global_short(particles);

	queue.enqueueNDRangeKernel(sum_kernel, cl::NullRange, global_short, cl::NullRange);
	queue.enqueueReadBuffer(buffer_x_1, CL_TRUE, 0, particles * sizeof(double), &x_1);
	queue.enqueueReadBuffer(buffer_y_1, CL_TRUE, 0, particles * sizeof(double), &y_1);

	for (int i = 0; i < particles; i++) {
		std::cout << "Expected sum x: " << expected_sum[i] << " Actual sum x: " << x_1[i] << "\n";
		std::cout << "Expected sum y: " << expected_sum[i] << " Actual sum y: " << y_1[i] << "\n";
	}
}