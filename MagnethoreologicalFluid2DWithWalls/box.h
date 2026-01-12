#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <fstream>
#include <string>
#include <math.h>
#include <cmath>

class Box {
private:
	int particles_;
	int length_;
	int dimensions_;
	std::vector<std::vector<double>> positions_;

	void InitialPositions();

	double CalculateInitialX(int position);

	double CalculateInitialY(int position);

	double CalculateInitialZ(int position);

	double GetX(int particle);
	double GetY(int particle);
	double GetZ(int particle);

public:
	Box(int particles, int length, int dimensions);
	~Box();

	void WritePositions(int iteration, double mason, double amplitude_relationship, int repetition, std::string tag);

	std::vector<double> ReturnX();

	std::vector<double> ReturnY();
	std::vector<double> ReturnZ();

	void SetX(double* x);

	void SetY(double* y);

	void SetZ(double* z);

	void ReadCsv(std::string path);
};