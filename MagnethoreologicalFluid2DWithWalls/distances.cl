#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void distances (global const double* X_0, global const double* Y_0, global const int* particle_0, global const int* particle_1, global const int* length, global const int* particles, global int* valid) {
	double r;
	double A = 2;
	double B = 10;
	double r_min = 1 - log(100.0)/B;
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

	r = distances[index];

	if (r < r_min){ *valid = 0; }
}