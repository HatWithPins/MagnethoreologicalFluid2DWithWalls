#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void validation (
	global int* valid, global double* x_0, global double* y_0, global double* x_1, global double* y_1, 
	global const double* original_delta_t, global double* delta_t, global double* t, global double* magnetic_field,
	global const double* mason, global const double* amplitude_relationship
) {
	int particle_0 = get_global_id(0);
	if (valid) {
		x_0[particle_0] = x_1[particle_0];
		y_0[particle_0] = y_1[particle_0];

		if (particle_0 == 0) {
			*t += *delta_t;
			*delta_t = *original_delta_t;

			magnetic_field[1] = 1/sqrt((*amplitude_relationship)*sin((*mason)*(*t))*(*amplitude_relationship)*sin((*mason)*(*t)) + 1);
			magnetic_field[0] = (*amplitude_relationship)*sin((*mason)*(*t))*magnetic_field[1];
		}
	}
	else if (particle_0 == 0) {
		*delta_t = *delta_t/2;
		*valid = 1;
	}
}