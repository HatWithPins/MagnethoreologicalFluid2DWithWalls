#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void sum (
	global double* x_0, global double* y_0, global double* z_0, global const int* dimensions, global const int* length, global const int* particles, global double* delta_t, 
	global double* forces_x, global double* forces_y, global double* forces_z, global const int* initial_indices_sum, global const int* last_indices_sum, global const int* matrix_size, 
	global int* valid, global double* x_1, global double* y_1, global double* z_1, global const int* mode, global const int* phase , global double* stress_array
	) {
		double r;
		double A = 1;
		double B = 100;
		double top_separation = 0.5;
		double top_repulsion = 0;
		double bottom_separation = 0.5;
		double bottom_repulsion = 0;
		double stress = 0;
		int particle_0 = get_global_id(0);
		int remainder_x;
		int initial_index_sum = initial_indices_sum[particle_0];
		int last_index_sum = last_indices_sum[particle_0];
		int index_sub = particle_0 - 1;
		double wall_velocity = (*length);

		x_1[particle_0] = 0;
		y_1[particle_0] = 0;
		z_1[particle_0] = 0;
		if (*mode == 1) {
			top_separation = (*length) - y_0[particle_0] + 0.5;
			top_repulsion = -A * exp(-B * (top_separation - 1));
			bottom_separation = y_0[particle_0] + 0.5;
			bottom_repulsion = A * exp(-B * (bottom_separation - 1));
			y_1[particle_0] += top_repulsion + bottom_repulsion;
		} else {
			top_separation = (*length) - z_0[particle_0] + 0.5;
			top_repulsion = -A * exp(-B * (top_separation - 1));
			bottom_separation = z_0[particle_0] + 0.5;
			bottom_repulsion = A * exp(-B * (bottom_separation - 1));
			z_1[particle_0] += top_repulsion + bottom_repulsion;

			if (((*mode == 2 || *mode == 3) && *phase == 1) || *phase == 2) {
				stress = z_0[particle_0] * wall_velocity / (*length);
			}
		}


		if (initial_index_sum < (*matrix_size)) {
			for (int i = initial_index_sum; i < last_index_sum; i++) {
				x_1[particle_0] += forces_x[i];
				y_1[particle_0] += forces_y[i];
				z_1[particle_0] += forces_z[i];
			}
		}
		if (particle_0 != 0) {
			for (int i = 0; i < particle_0; i++) {
				x_1[particle_0] -= forces_x[index_sub];
				y_1[particle_0] -= forces_y[index_sub];
				z_1[particle_0] -= forces_z[index_sub];

				index_sub += (*particles) - 2 - i;
			}
		}

		stress_array[particle_0] = x_1[particle_0]*z_0[particle_0];
		x_1[particle_0] += stress;
		x_1[particle_0] = (*delta_t)*x_1[particle_0];
		y_1[particle_0] = (*delta_t)*y_1[particle_0];
		z_1[particle_0] = (*delta_t)*z_1[particle_0];

		r = (sqrt(x_1[particle_0]*x_1[particle_0] + y_1[particle_0]*y_1[particle_0] + z_1[particle_0]*z_1[particle_0]));
		if (r > 0.03){ *valid = 0; }

		x_1[particle_0] += x_0[particle_0];
		y_1[particle_0] += y_0[particle_0];
		z_1[particle_0] += z_0[particle_0];

		remainder_x = x_1[particle_0]/(*length);
		x_1[particle_0] = x_1[particle_0] - remainder_x*(*length) + (*length)*(x_1[particle_0] < 0);
		if (*dimensions == 3) {
			int remainder_y = y_1[particle_0]/(*length);
			y_1[particle_0] = y_1[particle_0] - remainder_y*(*length) + (*length)*(y_1[particle_0] < 0);
		}
}