#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void distances (
	global const double* x_0, global const double* y_0, global const double* z_0, global const int* particle_0_array, global const int* particle_1_array, 
	global const int* length, global int* valid, global const int* dimensions
	) {
	double r = 1;
	double A = 2;
	double B = 10;
	double r_min = 1 - log(100.0)/B;
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
	} else {
		double distances[9];
		distances[0] = (sqrt(pown(x_0[particle_1] - x_0[particle_0], 2) + pown(y_0[particle_1] - y_0[particle_0], 2) + pown(z_0[particle_1]-z_0[particle_0],2)));
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
	}
	
	if (r < r_min){ *valid = 0; }
}