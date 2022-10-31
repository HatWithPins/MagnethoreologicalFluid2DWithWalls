#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void forces (
	global double* x_0, global double* y_0, global double* magnetic_field, global const int* length, global const int* particles, global const int* matrix_size,
	global const int* particle_0, global const int* particle_1, global double* forces_x, global double* forces_y
	) {
		double r;
		double A = 2;
		double B = 10;
		double x;
		double y;
		double dot_product;
		double distances[3];
		int idx = get_global_id(0);
		int particle_0 = particle_0[idx];
		int particle_1 = particle_1[idx];

		distances[0] = (sqrt(pown(X_0[particle_1] - X_0[particle_0],2) + pown(Y_0[particle_1] - Y_0[particle_0],2)));
		distances[1] = (sqrt(pown(X_0[particle_1] - X_0[particle_0] + (*length),2) + pown(Y_0[particle_1] - Y_0[particle_0],2)));
		distances[2] = (sqrt(pown(X_0[particle_1] - X_0[particle_0] - (*length),2) + pown(Y_0[particle_1] - Y_0[particle_0],2)));
		int index = 0;
		for(int j = 0; j < 9; j++) {
			if (distances[index] > distances[j]){ index = j; }
		}
		x = (*length)*((index == 1) - (index == 2)) + X_0[particula_1] - X_0[particula_0];
		y = Y_0[particula_1] - Y_0[particula_0];
		r = distances[index];
		x = x / r;
		y = y / r;

		dot_product = magnetic_field[0]*x + magnetic_field[1]*y;

		forces_x[idx] = ((1 - 5 * pown(dot_product, 2)) / (pown(r, 4)) + A * exp(-B * (r - 1))) * x + (2 * dot_product * field[0]) / (pown(r, 4));
		forces_y[idx] = ((1 - 5 * pown(dot_product, 2)) / (pown(r, 4)) + A * exp(-B * (r - 1))) * y + (2 * dot_product * field[1]) / (pown(r, 4));
}