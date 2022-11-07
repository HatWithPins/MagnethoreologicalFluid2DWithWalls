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
	int repetitions = 5;
	double delta_t = 0.001;

	//std::thread threads[5];
	Simulation(particles, box_length, 2.00, 1.00, 0.001, 0);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << "Total time: "
		<< duration.count() << " seconds" << endl;

	return 0;
}
