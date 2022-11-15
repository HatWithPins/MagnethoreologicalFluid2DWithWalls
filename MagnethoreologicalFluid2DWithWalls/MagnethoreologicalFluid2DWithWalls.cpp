// MagnethoreologicalFluid2DWithWalls.cpp : Defines the entry point for the application.
//

#include "MagnethoreologicalFluid2DWithWalls.h"

using namespace std;
using namespace std::chrono;

int main()
{
	auto start = high_resolution_clock::now();
	int particles = 400;
	int box_length = 44;
	double numbers[14] = { 0.001, 0.003, 0.004, 0.006, 0.009, 0.01, 0.04, 0.08, 0.1, 0.4, 0.6, 0.8, 1.0, 2.0 };
	double times[14] = { 6000, 6000, 6000, 6000, 3000, 3000, 3000, 3000, 1000, 1000, 1000, 1000, 1000, 1000 };
	int repetitions = 5;
	double delta_t = 0.001;

	std::thread threads[10];
	for (int i = 0; i < 14; i++) {
		for (int j = 0; j < 7; j++) {
			for (int k = 0; k < repetitions; k++) {
				threads[2*k] = std::thread(Simulation, particles, box_length, numbers[i], numbers[2*j], delta_t, k, times[i]);
				threads[2*k+1] = std::thread(Simulation, particles, box_length, numbers[i], numbers[2*j+1], delta_t, k, times[i]);
			}
			for (int k = 0; k < 2*repetitions; k++) {
				threads[k].join();
			}
		}
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << "Total time: "
		<< duration.count() << " seconds" << endl;

	return 0;
}
