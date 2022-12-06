#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void forces (
	global double* x_0, global double* y_0, global double* z_0, global const int* dimensions, global double* magnetic_field, global const int* length, global const int* particles, 
	global const int* matrix_size, global const int* particle_0_array, global const int* particle_1_array, global double* forces_x, global double* forces_y, global double* forces_z
	) {
		double r;
		double A = 1;
		double B = 100;
		double x;
		double y;
		double z;
		double dot_product;
		int idx = get_global_id(0);
		int particle_0 = particle_0_array[idx];
		int particle_1 = particle_1_array[idx];

		if (*dimensions == 2) {
			double distances[3];
			distances[0] = (sqrt(pown(x_0[particle_1] - x_0[particle_0],2) + pown(y_0[particle_1] - y_0[particle_0],2)));
			distances[1] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] + (*length),2) + pown(y_0[particle_1] - y_0[particle_0],2)));
			distances[2] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] - (*length),2) + pown(y_0[particle_1] - y_0[particle_0],2)));
			int index = 0;
			for(int j = 0; j < 3; j++) {
				if (distances[index] > distances[j]){ index = j; }
			}
			r = distances[index];

			x = (*length)*((index == 1) - (index == 2)) + x_0[particle_1] - x_0[particle_0];
			y = y_0[particle_1] - y_0[particle_0];
			x = x / r;
			y = y / r;

			dot_product = magnetic_field[0]*x + magnetic_field[1]*y;

			forces_x[idx] = ((1 - 5 * pown(dot_product, 2)) / (pown(r, 4)) + A * exp(-B * (r - 1))) * x + (2 * dot_product * magnetic_field[0]) / (pown(r, 4));
			forces_y[idx] = ((1 - 5 * pown(dot_product, 2)) / (pown(r, 4)) + A * exp(-B * (r - 1))) * y + (2 * dot_product * magnetic_field[1]) / (pown(r, 4));
			forces_z[idx] = 0;
		} else {
			double distances[9];
			distances[0] = (sqrt(pown(x_0[particle_1] - x_0[particle_0], 2) + pown(y_0[particle_1] - y_0[particle_0], 2) + pown(z_0[particle_1] - z_0[particle_0],2)));
			distances[1] = (sqrt(pown(x_0[particle_1] - x_0[particle_0], 2) + pown(y_0[particle_1] - y_0[particle_0] - (*length), 2) + pown(z_0[particle_1] - z_0[particle_0], 2)));
			distances[2] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] + (*length), 2) + pown(y_0[particle_1] - y_0[particle_0] - (*length), 2) + pown(z_0[particle_1] - z_0[particle_0], 2)));
			distances[3] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] + (*length), 2) + pown(y_0[particle_1] - y_0[particle_0], 2) + pown(z_0[particle_1] - z_0[particle_0], 2)));
			distances[4] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] + (*length), 2) + pown(y_0[particle_1] - y_0[particle_0] + (*length), 2) + pown(z_0[particle_1] - z_0[particle_0], 2)));
			distances[5] = (sqrt(pown(x_0[particle_1] - x_0[particle_0], 2) + pown(y_0[particle_1] - y_0[particle_0] + (*length), 2) + pown(z_0[particle_1] - z_0[particle_0], 2)));
			distances[6] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] - (*length), 2) + pown(y_0[particle_1] - y_0[particle_0] + (*length), 2) + pown(z_0[particle_1] - z_0[particle_0], 2)));
			distances[7] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] - (*length), 2) + pown(y_0[particle_1] - y_0[particle_0], 2) + pown(z_0[particle_1] - z_0[particle_0], 2)));
			distances[8] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] - (*length), 2) + pown(y_0[particle_1] - y_0[particle_0] - (*length), 2) + pown(z_0[particle_1] - z_0[particle_0], 2)));
			int index = 0;
			for(int j = 0; j < 9; j++) {
				if (distances[index] > distances[j]){ index = j; }
			}
			r = distances[index];

			x = (*length) * ((index == 2) + (index == 3) + (index == 4) - (index == 6) - (index == 7) - (index == 8)) + x_0[particle_1] - x_0[particle_0];
			y = (*length) * ((index == 4) + (index == 5) + (index == 6) - (index == 1) - (index == 2) - (index == 8)) + y_0[particle_1] - y_0[particle_0];
			z = z_0[particle_1] - z_0[particle_0];
			x = x / r;
			y = y / r;
			z = z / r;

			dot_product = magnetic_field[0]*x + magnetic_field[1]*y + magnetic_field[2]*z;

			forces_x[idx] = ((1 - 5 * pown(dot_product, 2)) / (pown(r, 4)) + A * exp(-B * (r - 1))) * x + (2 * dot_product * magnetic_field[0]) / (pown(r, 4));
			forces_y[idx] = ((1 - 5 * pown(dot_product, 2)) / (pown(r, 4)) + A * exp(-B * (r - 1))) * y + (2 * dot_product * magnetic_field[1]) / (pown(r, 4));
			forces_z[idx] = ((1 - 5 * pown(dot_product, 2)) / (pown(r, 4)) + A * exp(-B * (r - 1))) * z + (2 * dot_product * magnetic_field[2]) / (pown(r, 4));
		}
}