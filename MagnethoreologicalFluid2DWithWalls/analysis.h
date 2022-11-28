class Analysis {
private:
	double mason_;
	double amplitude_relationship_;
	int iteration_;
	int particles_;
	int length_;
	int dimensions_;
	int window_;
	double pi = 3.14159265359;
	double epsilon_ = 0.1;
	double micro_structure_separation = 1.1;
	double macro_structure_separation = 2.0;
	std::vector<std::vector<int>> micro_chains;
	std::vector<std::vector<int>> micro_sizes;
	std::vector<std::vector<double>> micro_linearities;
	std::vector<std::vector<double>> micro_means;
	std::vector<std::vector<int>> macro_chains;
	std::vector<std::vector<int>> macro_sizes;
	std::vector<std::vector<double>> macro_linearities;
	std::vector<std::vector<double>> macro_means;
	std::vector<double> times;

	int* Adjacency(double* x, double* y, double* z, double max_separation) {
		int* adjacency = new int[particles_ * particles_];
		double r;

		if (dimensions_ == 2) {
			double distances[3];
			for (int i = 0; i < particles_; i++) {
				for (int j = 0; j < particles_; j++) {
					distances[0] = (sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i], 2)));
					distances[1] = (sqrt(pow(x[j] - x[i] + (length_), 2) + pow(y[j] - y[i], 2)));
					distances[2] = (sqrt(pow(x[j] - x[i] - (length_), 2) + pow(y[j] - y[i], 2)));
					int index = 0;
					for (int k = 0; k < 3; k++) {
						if (distances[index] > distances[k]) { index = k; }
					}

					r = distances[index];

					adjacency[i * particles_ + j] = (r <= max_separation);
				}
			}
		}
		else if (dimensions_ == 3) {
			double distances[9];
			for (int i = 0; i < particles_; i++) {
				for (int j = 0; j < particles_; j++) {
					distances[0] = (sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i], 2) + pow(z[j]-z[i],2)));
					distances[1] = (sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i] - (length_), 2) + pow(z[j] - z[i], 2)));
					distances[2] = (sqrt(pow(x[j] - x[i] + (length_), 2) + pow(y[j] - y[i] - (length_), 2) + pow(z[j] - z[i], 2)));
					distances[3] = (sqrt(pow(x[j] - x[i] + (length_), 2) + pow(y[j] - y[i], 2) + pow(z[j] - z[i], 2)));
					distances[4] = (sqrt(pow(x[j] - x[i] + (length_), 2) + pow(y[j] - y[i] + (length_), 2) + pow(z[j] - z[i], 2)));
					distances[5] = (sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i] + (length_), 2) + pow(z[j] - z[i], 2)));
					distances[6] = (sqrt(pow(x[j] - x[i] - (length_), 2) + pow(y[j] - y[i] + (length_), 2) + pow(z[j] - z[i], 2)));
					distances[7] = (sqrt(pow(x[j] - x[i] - (length_), 2) + pow(y[j] - y[i], 2) + pow(z[j] - z[i], 2)));
					distances[8] = (sqrt(pow(x[j] - x[i] - (length_), 2) + pow(y[j] - y[i] - (length_), 2) + pow(z[j] - z[i], 2)));
					int index = 0;
					for (int k = 0; k < 9; k++) {
						if (distances[index] > distances[k]) { index = k; }
					}

					r = distances[index];

					adjacency[i * particles_ + j] = (r <= max_separation);
				}
			}
		}


		return adjacency;
	}

	std::vector<int> BFS(int* adjacency) {
		std::vector<int> chain;
		std::vector<int> list(1, 0);
		std::vector<int> visited_nodes(particles_, 0);

		for (int i = 0; i < particles_; i++) {
			chain.push_back(i);
		}

		bool visited = false;
		int node = 0;

		while (visited == false) {
			visited_nodes[list[0]] = 1;

			for (int i = 0; i < particles_; i++) {
				if (adjacency[particles_ * list[0] + i] == 1 && visited_nodes[i] == 0 && !(std::find(list.begin(), list.end(), i) != list.end())) {
					list.push_back(i);
				}
			}

			chain[list[0]] = min(node, chain[list[0]]);

			list.erase(list.begin());

			if (std::find(visited_nodes.begin(), visited_nodes.end(), 0) == visited_nodes.end()) {
				visited = true;
			}
			else {
				if (list.empty() == true) {
					for (int i = 0; i < particles_; i++) {
						if (visited_nodes[i] == 0) {
							list.push_back(i);
							node = i;
							break;
						}
					}
				}
			}
		}

		return chain;
	}

	std::vector<double> Linearity(double* x, double* y, double* z, std::vector<int> chain, std::vector<int> unique, std::vector<int> size) {
		std::vector<double> linearity;
		double xi, yi, Rx, Ry, Rz, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, lambda_1, lambda_2, lambda_3, I_max, I_min;
		int position = 0;
		double linear;

		for (int i : unique) {
			Rx = 0;
			Ry = 0;
			Rz = 0;
			Ixx = 0;
			Iyy = 0;
			Izz = 0;
			Ixy = 0;
			Ixz = 0;
			Iyz = 0;

			for (int j = 0; j < particles_; j++) {
				Rx += x[j] * (chain[j] == i);
				Ry += y[j] * (chain[j] == i);
				Rz += z[j] * (chain[j] == i);
			}

			Rx = Rx / size[position];
			Ry = Ry / size[position];
			Rz = Rz / size[position];

			for (int j = 0; j < particles_; j++) {
				if (dimensions_ == 3) {
					Ixx += (pow((Ry - y[j]), 2) + pow((Rz - z[j]), 2)) * (chain[j] == i);
					Iyy += (pow((Rx - x[j]), 2) + pow((Rz - z[j]), 2)) * (chain[j] == i);
					Izz += (pow((Rx - x[j]), 2) + pow((Ry - y[j]), 2)) * (chain[j] == i);
					Ixy -= (Ry - y[j]) * (Rx - x[j]) * (chain[j] == i);
					Ixz -= (Rz - z[j]) * (Rx - x[j]) * (chain[j] == i);
					Iyz -= (Ry - y[j]) * (Rz - z[j]) * (chain[j] == i);
				}
				else {
					Ixx += (pow((Ry - y[j]), 2)) * (chain[j] == i);
					Iyy += (pow((Rx - x[j]), 2)) * (chain[j] == i);
					Ixy -= (Ry - y[j]) * (Rx - x[j]) * (chain[j] == i);
					Ixz -= (Rz - z[j]) * (Rx - x[j]) * (chain[j] == i);
				}
			
			}

			if (dimensions_ == 3) {
				double p1 = pow(Ixy,2) + pow(Ixz,2) + pow(Iyz,2);
				if (p1 == 0) {
					lambda_1 = Ixx;
					lambda_2 = Iyy;
					lambda_3 = Izz;
				}
				else {
					double q = (Ixx + Iyy + Izz) / 3;
					double p2 = pow(Ixx - q, 2) + pow(Iyy - q, 2) + pow(Izz - q, 2) + 2 * p1;
					double p = sqrt(p2 / 6);
					double Bxx = (1 / p) * (Ixx - q);
					double Byy = (1 / p) * (Iyy - q);
					double Bzz = (1 / p) * (Izz - q);
					double Bxy = (1 / p) * Ixy;
					double Bxz = (1 / p) * Ixz;
					double Byz = (1 / p) * Iyz;
					double r = (Bxx*Byy*Bzz + 2*Bxy*Bxz*Byz - Bxx*Byz*Byz - Byy*Bxz*Bxz - Bzz*Bxy*Bxy) / 2;

					double phi;
					if (r <= -1) {
						phi = pi / 3;
					}
					else if (r >= 1) {
						phi = 0;
					}
					else {
						phi = acos(r) / 3;
					}

					lambda_1 = q + 2 * p * cos(phi);
					lambda_3 = q + 2 * p * cos(phi + (2 * pi / 3));
					lambda_2 = 3 * q - lambda_1 - lambda_3;
				}
			}
			else if (dimensions_ == 2) {
				lambda_1 = (Ixx + Iyy + sqrt((Ixx + Iyy) * (Ixx + Iyy) - 4 * (Ixx * Iyy - Ixy * Ixy))) / 2;
				lambda_2 = (Ixx + Iyy - sqrt((Ixx + Iyy) * (Ixx + Iyy) - 4 * (Ixx * Iyy - Ixy * Ixy))) / 2;
			}
			else {
				lambda_1 = 1.0;
				lambda_2 = 1.0;
			}
			

			I_max = (dimensions_ == 3) ? abs(max(lambda_1, lambda_2, lambda_3)) : abs(max(lambda_1, lambda_2));
			I_min = (dimensions_ == 3) ? abs(min(lambda_1, lambda_2, lambda_3)) : abs(min(lambda_1, lambda_2));

			linear = (sqrt(I_max) - sqrt(I_min)) / (sqrt(I_max) + sqrt(I_min));

			if (!isnan(linear)) {
				linearity.push_back(linear);
			}
			else {
				linearity.push_back(0);
			}

			position++;
		}

		return linearity;
	}

	void Averages() {
		int na;
		double average_na = 0;
		double sigma_na = 0;
		double average_size = 0;
		double sigma_size = 0;
		double average_linearity = 0;
		double sigma_linearity = 0;
		int index = iteration_ - 1;

		std::vector<double> average(6);
		average[1] = 0;
		average[2] = 0;
		average[5] = 0;


		na = micro_chains[index].size();

		for (int j = 0; j < na; j++) {
			if (micro_sizes[index][j] > 1) {
				average_na++;
				average_size += micro_sizes[index][j];
				average_linearity += micro_linearities[index][j];
			}
		}

		average[0] = average_na;
		average[2] = average_size / average_na;
		average[4] = average_linearity / average_na;

		for (int j = 0; j < na; j++) {
			if (micro_sizes[index][j] > 1) {
				sigma_na++;
				sigma_size += (micro_sizes[index][j] - average_size) * (micro_sizes[index][j] - average_size);
				sigma_linearity += (micro_linearities[index][j] - average_linearity) * (micro_linearities[index][j] - average_linearity);
			}
		}

		average[1] = 0;
		average[3] = sqrt(sigma_size / particles_);
		average[5] = sqrt(sigma_linearity / particles_);

		micro_means.push_back(average);

		//Macro structures averages
		average_na = 0;
		sigma_na = 0;
		average_size = 0;
		sigma_size = 0;
		average_linearity = 0;
		sigma_linearity = 0;

		average[1] = 0;
		average[2] = 0;
		average[5] = 0;


		na = macro_chains[index].size();

		for (int j = 0; j < na; j++) {
			if (macro_sizes[index][j] > 1) {
				average_na++;
				average_size += macro_sizes[index][j];
				average_linearity += macro_linearities[index][j];
			}
		}

		average[0] = average_na;
		average[2] = average_size / average_na;
		average[4] = average_linearity / average_na;

		for (int j = 0; j < na; j++) {
			if (macro_sizes[index][j] > 1) {
				sigma_na++;
				sigma_size += (macro_sizes[index][j] - average_size) * (macro_sizes[index][j] - average_size);
				sigma_linearity += (macro_linearities[index][j] - average_linearity) * (macro_linearities[index][j] - average_linearity);
			}
		}

		average[1] = 0;
		average[3] = sqrt(sigma_size / particles_);
		average[5] = sqrt(sigma_linearity / particles_);

		macro_means.push_back(average);
	}

	bool EndSimulation() {
		double micro_n = 0;
		double micro_average_na = 0;
		double micro_sigma_na = 0;
		double micro_average_size = 0;
		double micro_sigma_size = 0;
		double micro_average_linearity = 0;
		double micro_sigma_linearity = 0;

		for (int i = iteration_ - window_; i < iteration_ - 1; i++) {
			micro_n += micro_means[i][0];
			micro_average_size += micro_means[i][2] * micro_means[i][0];
			micro_average_linearity += micro_means[i][4] * micro_means[i][0];
		}

		micro_average_size = micro_average_size / micro_n;
		micro_average_linearity = micro_average_linearity / micro_n;
		micro_average_na = micro_n / iteration_;

		for (int i = iteration_ - window_; i < iteration_ - 1; i++) {
			micro_sigma_na += (micro_means[i][0] - micro_average_na) * (micro_means[i][0] - micro_average_na);
			micro_sigma_size += (micro_means[i][2] - micro_average_size) * (micro_means[i][2] - micro_average_size);
			micro_sigma_linearity += (micro_means[i][4] - micro_average_linearity) * (micro_means[i][4] - micro_average_linearity);
		}

		micro_sigma_na = sqrt(micro_sigma_na / micro_n);
		micro_sigma_size = sqrt(micro_sigma_size / micro_n);
		micro_sigma_linearity = sqrt(micro_sigma_linearity / micro_n);

		//Macro structures.
		double macro_n = 0;
		double macro_average_na = 0;
		double macro_sigma_na = 0;
		double macro_average_size = 0;
		double macro_sigma_size = 0;
		double macro_average_linearity = 0;
		double macro_sigma_linearity = 0;

		for (int i = iteration_ - window_; i < iteration_ - 1; i++) {
			macro_n += macro_means[i][0];
			macro_average_size += macro_means[i][2] * macro_means[i][0];
			macro_average_linearity += macro_means[i][4] * macro_means[i][0];
		}

		macro_average_size = macro_average_size / macro_n;
		macro_average_linearity = macro_average_linearity / macro_n;
		macro_average_na = macro_n / iteration_;

		for (int i = iteration_ - window_; i < iteration_ - 1; i++) {
			macro_sigma_na += (macro_means[i][0] - macro_average_na) * (macro_means[i][0] - macro_average_na);
			macro_sigma_size += (macro_means[i][2] - macro_average_size) * (macro_means[i][2] - macro_average_size);
			macro_sigma_linearity += (macro_means[i][4] - macro_average_linearity) * (macro_means[i][4] - macro_average_linearity);
		}

		macro_sigma_na = sqrt(macro_sigma_na / macro_n);
		macro_sigma_size = sqrt(macro_sigma_size / macro_n);
		macro_sigma_linearity = sqrt(macro_sigma_linearity / macro_n);

		return (micro_sigma_na / micro_average_na <= epsilon_) * (micro_sigma_size / micro_average_size <= epsilon_) * (micro_sigma_linearity / micro_average_linearity <= epsilon_)*
			(macro_sigma_na / macro_average_na <= epsilon_) * (macro_sigma_size / macro_average_size <= epsilon_) * (macro_sigma_linearity / macro_average_linearity <= epsilon_);
	}

public:

	Analysis(double mason, double amplitude_relationship, int particles, int length, int window, int dimensions) {
		mason_ = mason;
		amplitude_relationship_ = amplitude_relationship;
		particles_ = particles;
		length_ = length;
		window_ = window;
		dimensions_ = dimensions;
		iteration_ = 0;
	}
	~Analysis() {
		micro_chains.~vector();
		micro_linearities.~vector();
		micro_means.~vector();
		micro_sizes.~vector();
		macro_chains.~vector();
		macro_linearities.~vector();
		macro_means.~vector();
		macro_sizes.~vector();
	}

	bool PreAnalysis(double* x, double* y, double* z, double time) {
		iteration_++;
		int* micro_adjacency = Adjacency(x, y, z, micro_structure_separation);

		std::vector<int> micro_chain = BFS(micro_adjacency);

		std::vector<int> micro_unique(micro_chain.size());
		std::vector<int>::iterator it;

		it = std::unique_copy(micro_chain.begin(), micro_chain.end(), micro_unique.begin());
		std::sort(micro_unique.begin(), it);
		it = std::unique_copy(micro_unique.begin(), it, micro_unique.begin());
		micro_unique.resize(std::distance(micro_unique.begin(), it));

		std::vector<int> micro_length(micro_unique.size());

		for (size_t i = 0; i < micro_length.size(); ++i) {
			micro_length[i] = std::count(micro_chain.begin(), micro_chain.end(), micro_unique[i]);
		}

		std::vector<double> micro_linearity;
		micro_sizes.push_back(micro_length);

		micro_linearity = Linearity(x, y, z, micro_chain, micro_unique, micro_length);
		micro_linearities.push_back(micro_linearity);
		micro_chains.push_back(micro_unique);

		int* macro_adjacency = Adjacency(x, y, z, macro_structure_separation);

		std::vector<int> macro_chain = BFS(macro_adjacency);

		std::vector<int> macro_unique(macro_chain.size());

		it = std::unique_copy(macro_chain.begin(), macro_chain.end(), macro_unique.begin());
		std::sort(macro_unique.begin(), it);
		it = std::unique_copy(macro_unique.begin(), it, macro_unique.begin());
		macro_unique.resize(std::distance(macro_unique.begin(), it));

		std::vector<int> macro_length(macro_unique.size());

		for (size_t i = 0; i < macro_length.size(); ++i) {
			macro_length[i] = std::count(macro_chain.begin(), macro_chain.end(), macro_unique[i]);
		}

		std::vector<double> macro_linearity;
		macro_sizes.push_back(macro_length);

		macro_linearity = Linearity(x, y, z, macro_chain, macro_unique, macro_length);
		macro_linearities.push_back(macro_linearity);
		macro_chains.push_back(macro_unique);

		times.push_back(time);
		Averages();
		delete micro_adjacency;
		delete macro_adjacency;
		if (iteration_ >= window_ * 2) {
			return EndSimulation();
		}
		return false;
	}

	void WriteAnalysis(int repetition, std::string tag) {
		std::ofstream micro_analysis_file{ "analysis/micro_analysis-" + std::to_string(mason_) + "-" + std::to_string(amplitude_relationship_) + "-" + std::to_string(repetition) + "-" + tag + ".csv" };

		micro_analysis_file << "mason,amplitude,N,Na,sigma_Na,size,sigma_size,linearity,sigma_linearity\n";

		std::ofstream macro_analysis_file{ "analysis/macro_analysis-" + std::to_string(mason_) + "-" + std::to_string(amplitude_relationship_) + "-" + std::to_string(repetition) + "-" + tag + ".csv" };

		macro_analysis_file << "mason,amplitude,N,Na,sigma_Na,size,sigma_size,linearity,sigma_linearity\n";

		std::ofstream times_file{ "analysis/times-" + std::to_string(mason_) + "-" + std::to_string(amplitude_relationship_) + "-" + std::to_string(repetition) + "-" + tag + ".csv" };

		times_file << "t,micro_Na,micro_size,micro_linearity,macro_Na,macro_size,macro_linearity\n";

		double micro_n = 0;
		double micro_average_na = 0;
		double micro_sigma_na = 0;
		double micro_average_size = 0;
		double micro_sigma_size = 0;
		double micro_average_linearity = 0;
		double micro_sigma_linearity = 0;

		for (int i = 0; i < iteration_; i++) {
			micro_n += micro_means[i][0];
			micro_average_size += micro_means[i][2] * micro_means[i][0];
			micro_average_linearity += micro_means[i][4] * micro_means[i][0];
		}

		micro_average_size = micro_average_size / micro_n;
		micro_average_linearity = micro_average_linearity / micro_n;
		micro_average_na = micro_n / iteration_;

		for (int i = 0; i < iteration_; i++) {
			micro_sigma_na += (micro_means[i][0] - micro_average_na) * (micro_means[i][0] - micro_average_na);
			micro_sigma_size += (micro_means[i][2] - micro_average_size) * (micro_means[i][2] - micro_average_size);
			micro_sigma_linearity += (micro_means[i][4] - micro_average_linearity) * (micro_means[i][4] - micro_average_linearity);
		}

		micro_sigma_na = sqrt(micro_sigma_na / micro_n);
		micro_sigma_size = sqrt(micro_sigma_size / micro_n);
		micro_sigma_linearity = sqrt(micro_sigma_linearity / micro_n);

		//Macro structures.
		double macro_n = 0;
		double macro_average_na = 0;
		double macro_sigma_na = 0;
		double macro_average_size = 0;
		double macro_sigma_size = 0;
		double macro_average_linearity = 0;
		double macro_sigma_linearity = 0;

		for (int i = 0; i < iteration_; i++) {
			macro_n += macro_means[i][0];
			macro_average_size += macro_means[i][2] * macro_means[i][0];
			macro_average_linearity += macro_means[i][4] * macro_means[i][0];
		}

		macro_average_size = macro_average_size / macro_n;
		macro_average_linearity = macro_average_linearity / macro_n;
		macro_average_na = macro_n / iteration_;

		for (int i = 0; i < iteration_; i++) {
			macro_sigma_na += (macro_means[i][0] - macro_average_na) * (macro_means[i][0] - macro_average_na);
			macro_sigma_size += (macro_means[i][2] - macro_average_size) * (macro_means[i][2] - macro_average_size);
			macro_sigma_linearity += (macro_means[i][4] - macro_average_linearity) * (macro_means[i][4] - macro_average_linearity);
		}

		macro_sigma_na = sqrt(macro_sigma_na / macro_n);
		macro_sigma_size = sqrt(macro_sigma_size / macro_n);
		macro_sigma_linearity = sqrt(macro_sigma_linearity / macro_n);

		for (int i = 0; i < iteration_; i++) {
			times_file << times[i] << "," << micro_means[i][0] << "," << micro_means[i][2] << "," << micro_means[i][4] << "," 
				<< macro_means[i][0] << "," << macro_means[i][2] << "," << macro_means[i][4] << "\n";
		}

		micro_analysis_file << mason_ << "," << amplitude_relationship_ << "," << micro_n << "," 
			<< micro_average_na << "," << micro_sigma_na << "," << micro_average_size << "," << micro_sigma_size << "," << micro_average_linearity << "," << micro_sigma_linearity 
			<< "\n";

		macro_analysis_file << mason_ << "," << amplitude_relationship_ << "," << macro_n << ","
			<< macro_average_na << "," << macro_sigma_na << "," << macro_average_size << "," << macro_sigma_size << "," << macro_average_linearity << "," << macro_sigma_linearity
			<< "\n";

		times_file.close();
		micro_analysis_file.close();
		macro_analysis_file.close();
	}
};