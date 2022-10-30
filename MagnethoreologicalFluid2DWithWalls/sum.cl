#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void sum (global const double* x_0, global const double* y_0, global const int* length, global const int* particles, global const double* delta_t, global double* forces_x, global double forces_y global int* valid, global double* x_1, global double* y_1) {
	double r;
	int particle_0 = get_global_id(0);
	int remainder_x;
	int matrix_size = (*particles)*(*particles - 1)/2;
	int initial_index_sum = 0;
	int last_index_sum = 0;
	int index_sub = particle_0 - 1;

	if (particle_0 != 0) {
		for (int i = 1; i <= particle_0; i++) {
			initial_index_sum += (*particles) - i;
		}
	}
	for (int i = 1; i <= particle_0 + 1; i++) {
		last_index_sum += (*particles) - i;
	}

	x_1[particle_0] = 0;
	y_1[particle_0] = 0;

	if (initial_index_sum < matrix_size) {
		for (int i = initial_index_sum; i < last_index_sum; i++) {
			x_1[particle_0] += forces_x[i];
			y_1[particle_0] += forces_y[i];
		}
	}
	if (particle_0 != 0) {
		for (int i = 0; i < particle_0; i++) {
			x_1[particle_0] -= forces_x[index_sub];
			y_1[particle_0] -= forces_y[index_sub];

			index_sub += (*particles) - i;
		}
	}

	x_1[particle_0] = (*delta_t)*x_1[particle_0];
	y_1[particle_0] = (*delta_t)*y_1[particle_0];

	r = (sqrt(x_1[particle_0]*x_1[particle_0] + y_1[particle_0]*y_1[particle_0]);
	valid[particle_0] = (r <= 0.03);

	x_1[particle_0] += x_0[particle_0];
	y_1[particle_0] += y_0[particle_0];

	remainder_x = x_1[particle_0]/(*length);
	x_1[particle_0] = x_1[particle_0] - remainder_x*(*length) + (*length)*(x_1[particle_0] < 0)
}