// MagnethoreologicalFluid2DWithWalls.cpp : Defines the entry point for the application.
//

#include "MagnethoreologicalFluid2DWithWalls.h"

using namespace std;
using namespace std::chrono;

int Length(int particles, int dimensions, double concentration) {
	double pi = 3.14159265359;
	int length = 0;

	if (dimensions == 3) {
		length = cbrt(particles * (pi / 6) / concentration);
	}
	else if (dimensions == 2) {
		double concentration2d = (pi / 4) * pow(cbrt(6*concentration/pi), 2);
		length = sqrt(particles / concentration2d);
	}
	
	return length;
}

void SimulationThread(int repetition, int particles, double mason, double ar, int dimensions, double concentration, double field_direction, int keep_positions, int load_positions, double creep_time) {
	
	double times[3];
	int phases;
	if (mason > 0.0) {
		phases = 3;
		times[0] = 500.0;
		times[1] = 1000.0;
		times[2] = 1050.0;
	}
	else {
		phases = 2;
		times[0] = 500.0;
		times[1] = 550.0;
		times[2] = 550.0;
	}

	int boxLength = Length(particles, dimensions, concentration);
	double deltaT = 0.001;

	Simulation(field_direction, phases, particles, dimensions, boxLength, mason, ar, deltaT, repetition, times, keep_positions, load_positions, creep_time);
}

int main(int argc, char** argv)
{
	auto start = high_resolution_clock::now();

	string argument;
	string path;
	size_t pos;
	size_t check;
	int repetitions = 1;
	int particles = 400;
	double ma = 0.0;
	double ar = 0.0;
	double concentration = 0.07;
	double creep_time = 0.0;
	int dimensions = 2;
	//Angle of the perturbation with respect of x axis. 
	double field_direction = 0;
	int keep_positions = 0;
	int load_positions = 0;
	int expectedArguments = 11;
	vector<string> expectedArgumentsList = { 
		"repetitions=", 
		"particles=",
		"concentration=",
		"ma=",
		"ar=", 
		"dimensions=",
		"field_direction=",
		"keep_positions=",
		"load_positions=",
		"creep_time="
	};

	if (argc > expectedArguments)
	{
		cout << "Error, expected " << expectedArguments - 1 << ", but received " << argc - 1 << endl;
		return -1;
	}

	try
	{
		for (int i = 1; i < argc; i++)
		{
			argument = argv[i];
			check = argument.find(expectedArgumentsList[i - 1]);
			if (check < 0 || check > argument.size())
			{
				vector<string> exceptionVector = { argument, expectedArgumentsList[i - 1] };
				throw(exceptionVector);
			}
			pos = argument.find("=");

			try
			{
				if (argument.substr(0, pos) == "repetitions")
				{
					repetitions = stoi(argument.substr(pos + 1));
					if (repetitions < 0 || repetitions > 5) {
						cout << "Error, repetitions must be between 0 and 5, but received " << particles << endl;
						return -1;
					}
				}
				else if (argument.substr(0, pos) == "particles")
				{
					particles = stoi(argument.substr(pos + 1));
					if (particles < 2) {
						cout << "Error, particles must be greater than 1, but received " << particles << endl;
						return -1;
					}
				}
				else if (argument.substr(0, pos) == "concentration")
				{
					concentration = stod(argument.substr(pos + 1));
					if (concentration <= 0) {
						cout << "Error, concentration must be greater than 0, but received " << concentration << endl;
						return -1;
					}
				}
				else if (argument.substr(0, pos) == "ma")
				{
					ma = stod(argument.substr(pos + 1));
					if (ma < 0) {
						cout << "Error, ma must be greater or equal than 0, but received " << ma << endl;
						return -1;
					}
				}
				else if (argument.substr(0, pos) == "ar")
				{
					ar = stod(argument.substr(pos + 1));
					if (ar < 0) {
						cout << "Error, ar must be greater than 0, but received " << ar << endl;
						return -1;
					}
				}
				else if (argument.substr(0, pos) == "dimensions")
				{
					dimensions = stoi(argument.substr(pos + 1));
					if (dimensions > 3 || dimensions < 2) {
						cout << "Error, dimensions must be either 2 or 3, but received " << dimensions << endl;
						return -1;
					}
				}
				else if (argument.substr(0, pos) == "field_direction")
				{
					field_direction = stod(argument.substr(pos + 1));
					if (field_direction > 90 || field_direction < 0) {
						cout << "Error, field_direction must be either 0 or 90, but received " << field_direction << endl;
						return -1;
					}
				}
				else if (argument.substr(0, pos) == "keep_positions")
				{
					keep_positions = stoi(argument.substr(pos + 1));
					if (keep_positions > 1 || keep_positions < 0) {
						cout << "Error, keep_positions must be either 0 or 1, but received " << keep_positions << endl;
						return -1;
					}
				}
				else if (argument.substr(0, pos) == "load_positions")
				{
					load_positions = stoi(argument.substr(pos + 1));
					if (load_positions > 1 || load_positions < 0) {
						cout << "Error, load_positions must be either 0 or 1, but received " << load_positions << endl;
						return -1;
					}
				}
				else if (argument.substr(0, pos) == "creep_time")
				{
					creep_time = stod(argument.substr(pos + 1));
					if (ar < 0) {
						cout << "Error, creep_time must be greater than 0, but received " << creep_time << endl;
						return -1;
					}
				}
			}
			catch (const std::exception& e)
			{
				std::cerr << e.what() << '\n';
				return -1;
			}
		}
	}
	catch (vector<string> errorVector)
	{
		cout << "Error, expected " << errorVector[1] << "something, but received " << errorVector[0] << endl;
		return -1;
	}


	for (int i = 0; i < repetitions; i++) {
		SimulationThread( i, particles, ma, ar, dimensions, concentration, field_direction, keep_positions, load_positions, creep_time);
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	std::cout << "Total time: "
		<< duration.count() << " seconds" << endl;

	return 0;
}
