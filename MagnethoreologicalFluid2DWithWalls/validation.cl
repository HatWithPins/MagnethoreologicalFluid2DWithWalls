#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void validation (
	global int* valid, global double* x_0, global double* y_0, global double* x_1, global double* y_1, 
	global const double* original_delta_t, global double* delta_t global double* t, global double* magnetic_field,
	global const double* mason, global const double* amplitude_relationship
) {
	if (valid) {
		t += delta_t;
		delta_t = original_delta_t;
	}
	else {
	}
}