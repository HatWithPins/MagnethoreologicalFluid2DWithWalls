class Box {
private:
	int m_particles;
	int m_length;
	std::vector<std::vector<double>> positions;

	void InitialPositions() {
		std::vector<int> initial_positions{};
		std::vector<double> x{};
		std::vector<double> y{};

		for (int i = 0; i < m_length * m_length; i++) {
			initial_positions.push_back(i);
		}

		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::shuffle(initial_positions.begin(), initial_positions.end(), std::default_random_engine(seed));

		for (int i = 0; i < m_particles; i++) {
			x.push_back(CalculateInitialX(initial_positions[i]));
			y.push_back(CalculateInitialY(initial_positions[i]));
		}
		positions.push_back(x);
		positions.push_back(y);
	}

	double CalculateInitialX(int position) {
		return position % m_length + 0.5;
	}

	double CalculateInitialY(int position) {
		int division = position / m_length;
		return division + 0.5;
	}

	double GetX(int particle) {
		return positions[0][particle];
	}
	double GetY(int particle) {
		return positions[1][particle];
	}

public:
	Box(int particles, int length) {
		m_particles = particles;
		m_length = length;
		InitialPositions();
	}

	void WritePositions(int iteration, double mason, int repetition, std::string disturbed, double amplitude = 0.0) {
		std::ofstream file{ "positions/positions-" + std::to_string(mason) + "-" + std::to_string(amplitude) + "-" + std::to_string(repetition) + "-" + std::to_string(iteration) + "-" + disturbed + ".csv" };

		file << "x,y\n";

		for (int i = 0; i < m_particles; i++) {
			file << GetX(i) << "," << GetY(i) << "\n";
		}
		file.close();
	}

	std::vector<double> ReturnX() {
		return positions[0];
	}

	std::vector<double> ReturnY() {
		return positions[1];
	}

	void SetX(double* x) {
		for (int i = 0; i < m_particles; i++) {
			positions[0][i] = x[i];
		}
	}

	void SetY(double* y) {
		for (int i = 0; i < m_particles; i++) {
			positions[1][i] = y[i];
		}
	}
};