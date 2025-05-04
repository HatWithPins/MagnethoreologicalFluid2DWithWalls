#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void validation (
	global int* valid, global double* x_0, global double* y_0, global double* z_0, global double* x_1, global double* y_1, global double* z_1,
	global const double* original_delta_t, global double* delta_t, global double* t, global double* magnetic_field,
	global const double* mason, global const double* amplitude_relationship, global const int* dimensions, global const int* mode, global const int* phase,
	global double* stress_array, global double* stress, global const int* particles, global const int* length, global const int* field_direction
) {
	int particle_0 = get_global_id(0);
	if (*valid == 1) {
		x_0[particle_0] = x_1[particle_0];
		y_0[particle_0] = y_1[particle_0];
		z_0[particle_0] = z_1[particle_0];
		double volume = (*length)*(*length)*(*length);

		if (particle_0 == 0) {
			*t += *delta_t;
			*delta_t = *original_delta_t;

			if (*dimensions == 2) {
				magnetic_field[2] = 0;
				magnetic_field[1] = 1/sqrt((*amplitude_relationship)*sin((*mason)*(*t))*(*amplitude_relationship)*sin((*mason)*(*t)) + 1);
				magnetic_field[0] = (*amplitude_relationship)*sin((*mason)*(*t))*magnetic_field[1];
			} else if (*field_direction == 1) {
				magnetic_field[0] = 0;
				magnetic_field[2] = 1/sqrt((*amplitude_relationship)*sin((*mason)*(*t))*(*amplitude_relationship)*sin((*mason)*(*t)) + 1);
				magnetic_field[1] = (*amplitude_relationship)*sin((*mason)*(*t))*magnetic_field[2];

				if (*mode == 1) {
					*stress = 0;
					for (int i = 0; i < *particles; i++) {
						*stress += stress_array[i]; 
					}
					*stress = -(*stress) / volume;
				}
			} else if (*field_direction == 0) {
				magnetic_field[1] = 0;
				magnetic_field[2] = 1/sqrt((*amplitude_relationship)*sin((*mason)*(*t))*(*amplitude_relationship)*sin((*mason)*(*t)) + 1);
				magnetic_field[0] = (*amplitude_relationship)*sin((*mason)*(*t))*magnetic_field[2];

				if (*mode == 1) {
					*stress = 0;
					for (int i = 0; i < *particles; i++) {
						*stress += stress_array[i]; 
					}
					*stress = -(*stress) / volume;
				}
			} 
		}
	}
	else if (particle_0 == 0) {
		*delta_t = *delta_t/2;
		*valid = 1;
	}
}