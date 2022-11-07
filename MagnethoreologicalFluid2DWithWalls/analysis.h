class Analysis {
private:
	double mason_;
	double amplitude_relationship_;
	int iteration_;
	int particles_;
	int length_;
	int window_;
	double epsilon_ = 0.1;
	double micro_structure_separation = 1.1;
	double macro_structure_separation = 1.5;
	std::vector<std::vector<int>> chains;
	std::vector<std::vector<int>> sizes;
	std::vector<std::vector<double>> linearities;
	std::vector<std::vector<double>> means;
	std::vector<std::vector<double>> first_moving_average;
	std::vector<std::vector<double>> second_moving_average;
	std::vector<double> times;

	int* Adjacency(double* x, double* y, double max_separation) {
		int* adjacency = new int[particles_ * particles_];
		double r;
		double distances[3];

		for (int i = 0; i < particles_; i++) {
			for (int j = 0; j < particles_; j++) {
				distances[0] = (sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i], 2)));
				distances[1] = (sqrt(pow(x[j] - x[i] + (length_), 2) + pow(y[j] - y[i], 2)));
				distances[2] = (sqrt(pow(x[j] - x[i] - (length_), 2) + pow(y[j] - y[i], 2)));
				int index = 0;
				for (int k = 0; k < 9; k++) {
					if (distances[index] > distances[k]) { index = k; }
				}

				r = distances[index];

				adjacency[i * particles_ + j] = (r <= max_separation);
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

	std::vector<double> Linearity(double* x, double* y, std::vector<int> chain, std::vector<int> unique, std::vector<int> size) {
		std::vector<double> linearity;
		double xi, yi, Rx, Ry, Ixx, Iyy, Ixy, lambda_1, lambda_2, I_max, I_min;
		int position = 0;
		double linear;

		for (int i : unique) {
			Rx = 0;
			Ry = 0;
			Ixx = 0;
			Iyy = 0;
			Ixy = 0;

			for (int j = 0; j < particles_; j++) {
				Rx += x[j] * (chain[j] == i);
				Ry += y[j] * (chain[j] == i);
			}

			Rx = Rx / size[position];
			Ry = Ry / size[position];

			for (int j = 0; j < particles_; j++) {
				Ixx += (Ry - y[j]) * (Ry - y[j]) * (chain[j] == i);
				Iyy += (Rx - x[j]) * (Rx - x[j]) * (chain[j] == i);
				Ixy -= (Ry - y[j]) * (Rx - x[j]) * (chain[j] == i);
			}

			lambda_1 = (Ixx + Iyy + sqrt((Ixx + Iyy) * (Ixx + Iyy) - 4 * (Ixx * Iyy - Ixy * Ixy))) / 2;
			lambda_2 = (Ixx + Iyy - sqrt((Ixx + Iyy) * (Ixx + Iyy) - 4 * (Ixx * Iyy - Ixy * Ixy))) / 2;

			I_max = max(lambda_1, lambda_2);
			I_min = min(lambda_1, lambda_2);

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


		na = chains[index].size();

		for (int j = 0; j < na; j++) {
			if (sizes[index][j] > 1) {
				average_na++;
				average_size += sizes[index][j];
				average_linearity += linearities[index][j];
			}
		}

		average[0] = average_na;
		average[2] = average_size / average_na;
		average[4] = average_linearity / average_na;

		for (int j = 0; j < na; j++) {
			if (sizes[index][j] > 1) {
				sigma_na++;
				sigma_size += (sizes[index][j] - average_size) * (sizes[index][j] - average_size);
				sigma_linearity += (linearities[index][j] - average_linearity) * (linearities[index][j] - average_linearity);
			}
		}

		average[1] = 0;
		average[3] = sqrt(sigma_size / particles_);
		average[5] = sqrt(sigma_linearity / particles_);

		means.push_back(average);
	}

	bool EndSimulation() {
		double n = 0;
		double average_na = 0;
		double sigma_na = 0;
		double average_size = 0;
		double sigma_size = 0;
		double average_linearity = 0;
		double sigma_linearity = 0;

		for (int i = iteration_ - window_; i < iteration_ - 1; i++) {
			n += means[i][0];
			average_size += means[i][2] * means[i][0];
			average_linearity += means[i][4] * means[i][0];
		}

		average_size = average_size / n;
		average_linearity = average_linearity / n;
		average_na = n / iteration_;

		for (int i = iteration_ - window_; i < iteration_ - 1; i++) {
			sigma_na += (means[i][0] - average_na) * (means[i][0] - average_na);
			sigma_size += (means[i][2] - average_size) * (means[i][2] - average_size);
			sigma_linearity += (means[i][4] - average_linearity) * (means[i][4] - average_linearity);
		}

		sigma_na = sqrt(sigma_na / n);
		sigma_size = sqrt(sigma_size / n);
		sigma_linearity = sqrt(sigma_linearity / n);

		return (sigma_na / average_na <= epsilon_) * (sigma_size / average_size <= epsilon_) * (sigma_linearity / average_linearity <= epsilon_);
	}

public:

	Analysis(double mason, double amplitude_relationship, int particles, int length, int window) {
		mason_ = mason;
		amplitude_relationship_ = amplitude_relationship;
		length_ = length;
		window_ = window;
		iteration_ = 0;
	}

	bool PreAnalysis(double* x, double* y, double time) {
		iteration_++;
		int* adjacency = new int[particles_ * particles_];
		adjacency = Adjacency(x, y, micro_structure_separation);

		std::vector<int> chain = BFS(adjacency);

		std::vector<int> unique(chain.size());
		std::vector<int>::iterator it;

		it = std::unique_copy(chain.begin(), chain.end(), unique.begin());
		std::sort(unique.begin(), it);
		it = std::unique_copy(unique.begin(), it, unique.begin());
		unique.resize(std::distance(unique.begin(), it));

		std::vector<int> length(unique.size());

		for (size_t i = 0; i < length.size(); ++i) {
			length[i] = std::count(chain.begin(), chain.end(), unique[i]);
		}

		std::vector<double> linearity;
		sizes.push_back(length);

		linearity = Linearity(x, y, chain, unique, length);
		linearities.push_back(linearity);
		chains.push_back(unique);

		times.push_back(time);
		Averages();
		if (iteration_ >= window_ * 2) {
			return EndSimulation();
		}
		return false;
	}

	void WriteAnalysis(int repetition, std::string tag) {
		std::ofstream analysis_file{ "analysis/analysis-" + std::to_string(mason_) + "-" + std::to_string(amplitude_relationship_) + "-" + std::to_string(repetition) + "-" + tag + ".csv" };

		analysis_file << "mason,amplitude,N,Na,sigma_Na,size,sigma_size,linearity,sigma_linearity\n";

		std::ofstream times_file{ "analysis/times-" + std::to_string(mason_) + "-" + std::to_string(amplitude_relationship_) + "-" + std::to_string(repetition) + "-" + tag + ".csv" };

		times_file << "t,Na,size,linearity\n";

		double n = 0;
		double average_na = 0;
		double sigma_na = 0;
		double average_size = 0;
		double sigma_size = 0;
		double average_linearity = 0;
		double sigma_linearity = 0;

		for (int i = iteration_ - window_; i < iteration_ - 1; i++) {
			n += means[i][0];
			average_size += means[i][2] * means[i][0];
			average_linearity += means[i][4] * means[i][0];
			times_file << times[i] << "," << means[i][0] << "," << means[i][2] << "," << means[i][4] << "\n";
		}

		average_size = average_size / n;
		average_linearity = average_linearity / n;
		average_na = n / iteration_;

		for (int i = iteration_ - window_; i < iteration_ - 1; i++) {
			sigma_na += (means[i][0] - average_na) * (means[i][0] - average_na);
			sigma_size += (means[i][2] - average_size) * (means[i][2] - average_size);
			sigma_linearity += (means[i][4] - average_linearity) * (means[i][4] - average_linearity);
		}

		sigma_na = sqrt(sigma_na / n);
		sigma_size = sqrt(sigma_size / n);
		sigma_linearity = sqrt(sigma_linearity / n);

		analysis_file << mason_ << "," << amplitude_relationship_ << "," << n << "," << average_na << "," << sigma_na << "," << average_size << "," << sigma_size << "," << average_linearity << "," << sigma_linearity << "\n";

		times_file.close();
		analysis_file.close();
	}
};