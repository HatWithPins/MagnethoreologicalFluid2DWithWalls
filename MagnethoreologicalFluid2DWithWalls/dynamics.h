#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

void Simulation(double fieldDirection, int phases, int particles, int dimensions, int length, double mason, double amplitude_relationship, double original_delta_t, int repetition, double max_times[3], bool keep_positions=false, bool load_positions=false, double creep_time=0.0);
#endif