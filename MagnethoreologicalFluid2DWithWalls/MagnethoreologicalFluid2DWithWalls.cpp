﻿// MagnethoreologicalFluid2DWithWalls.cpp : Defines the entry point for the application.
//

#include "MagnethoreologicalFluid2DWithWalls.h"

using namespace std;
using namespace std::chrono;

int Length(int particles, int dimensions, double concentration) {
	double pi = 3.14159265359;
	int length = 0;

	if (dimensions == 3) {
		length = cbrt(particles / concentration);
	}
	else if (dimensions == 2) {
		double concentration_2d = (pi / 4) * pow(cbrt(6*concentration/pi), 2);
		length = sqrt(particles / concentration_2d);
	}
	
	return length;
}

int main()
{
	auto start = high_resolution_clock::now();

	std::thread threads[6];
	int repetitions = 5;
	double numbers[6] = { 0.4, 0.6, 0.8, 1.0, 1.5, 2.0 };

	for (int rep; rep < repetitions; rep++) {
		for (int i = 0; i < 5; i++) { //Mode
			for (int j = 0; j < 6; j++) { //RA
				for (int k = 0; k < 3; k++) { //Ma

					int mode = i + 2;
					int particles = 400;
					double concentration = 0.07;
					double times[3];
					int phases;
					int dimensions;

				switch (mode) {
				case 1:
					dimensions = 2;
					times[0] = 500;
					times[1] = 500;
					times[2] = 500;
					phases = 1;
					break;
				case 2:
					dimensions = 3;
					times[0] = 500;
					times[1] = 550;
					times[2] = 550;
					phases = 2;
					break;
				case 3:
					dimensions = 3;
					times[0] = 500;
					times[1] = 550;
					times[2] = 550;
					phases = 2;
					break;
				case 4:
					dimensions = 3;
					times[0] = 500;
					times[1] = 1000;
					times[2] = 1050;
					phases = 3;
					break;
				case 5:
					dimensions = 3;
					times[0] = 500;
					times[1] = 1000;
					times[2] = 1050;
					phases = 3;
					break;
				case 6:
					dimensions = 3;
					times[0] = 500;
					times[1] = 1000;
					times[2] = 1050;
					phases = 3;
					break;
				}

					int box_length = Length(particles, dimensions, concentration);
					double delta_t = 0.001;

					threads[2 * k] = std::thread(Simulation, mode, phases, particles, dimensions, box_length, numbers[2 * k], numbers[j], delta_t, rep, times);
					threads[2 * k + 1] = std::thread(Simulation, mode, phases, particles, dimensions, box_length, numbers[2 * k + 1], numbers[j], delta_t, rep, times);
				}
				for (int l = 0; l < 6; l++) {
					threads[l].join();
				}
			}
		}
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	std::cout << "Total time: "
		<< duration.count() << " seconds" << endl;

	return 0;
}
