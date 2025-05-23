#include <vector>
#include <string>

class Analysis {
private:
	double mason_;
	double amplitude_relationship_;
	int iteration_;
	int particles_;
	int length_;
	int dimensions_;
	double field_direction_;
	int window_;
	double pi = 3.14159265359;
	double epsilon_ = 0.1;
	double micro_structure_separation = 1.1;
	double macro_structure_separation = 2.0;
	double magnetic_field[3];
	std::vector<std::vector<int>> micro_chains;
	std::vector<std::vector<int>> micro_sizes;
	std::vector<std::vector<double>> micro_linearities;
	std::vector<std::vector<double>> micro_means;
	std::vector<std::vector<int>> macro_chains;
	std::vector<std::vector<int>> macro_sizes;
	std::vector<std::vector<double>> macro_linearities;
	std::vector<std::vector<double>> macro_means;
	std::vector<double> times;
	std::vector<double> stress_times;
	std::vector<double> stress_vector;
	std::vector<std::vector<double>> connectivity_vector;

	int* Adjacency(double* x, double* y, double* z, double max_separation);
	std::vector<int> BFS(int* adjacency);
	std::vector<double> Linearity(double* x, double* y, double* z, std::vector<int> chain, std::vector<int> unique, std::vector<int> size);
	void Averages();
	bool EndSimulation();

public:
	Analysis(double mason, double amplitude_relationship, int particles, int length, int window, int dimensions, double field_direction);
	~Analysis();
	bool PreAnalysis(double* x, double* y, double* z, double time);
	void Connectivity(double* x, double* y, double* z);
	void WriteConnectivity(int repetition, std::string tag);
	void RecordStress(double t, double stress);
	void WriteStress(int repetition, std::string tag);
	void WriteAnalysis(int repetition, std::string tag);
};