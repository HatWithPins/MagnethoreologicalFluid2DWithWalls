#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void distances (
	global const double* x_0, global const double* y_0, global const int* particle_0_array, global const int* particle_1_array, 
	global const int* length, global int* valid
	) {
	double r;
	double A = 2;
	double B = 10;
	double r_min = 1 - log(100.0)/B;
	double distances[3];
	int idx = get_global_id(0);
	int particle_0 = particle_0_array[idx];
	int particle_1 = particle_1_array[idx];

	distances[0] = (sqrt(pown(x_0[particle_1] - x_0[particle_0],2) + pown(y_0[particle_1] - y_0[particle_0],2)));
	distances[1] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] + (*length),2) + pown(y_0[particle_1] - y_0[particle_0],2)));
	distances[2] = (sqrt(pown(x_0[particle_1] - x_0[particle_0] - (*length),2) + pown(y_0[particle_1] - y_0[particle_0],2)));
	int index = 0;
	for(int j = 0; j < 3; j++) {
		if (distances[index] > distances[j]){ index = j; }
	}

	r = distances[index];

	if (r < r_min){ *valid = 0; }
}