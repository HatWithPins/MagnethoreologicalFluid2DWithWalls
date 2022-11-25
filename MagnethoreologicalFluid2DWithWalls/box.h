class Box {
private:
	int particles_;
	int length_;
	int dimensions_;
	std::vector<std::vector<double>> positions_;

	void InitialPositions() {
		std::vector<int> initial_positions{};
		std::vector<double> x{};
		std::vector<double> y{};
		std::vector<double> z{};

		for (int i = 0; i < pow(length_, dimensions_); i++) {
			initial_positions.push_back(i);
		}

		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::shuffle(initial_positions.begin(), initial_positions.end(), std::default_random_engine(seed));

		for (int i = 0; i < particles_; i++) {
			x.push_back(CalculateInitialX(initial_positions[i]));
			y.push_back(CalculateInitialY(initial_positions[i]));
			z.push_back(CalculateInitialZ(initial_positions[i]));
		}
		positions_.push_back(x);
		positions_.push_back(y);
		positions_.push_back(z);
	}

	double CalculateInitialX(int position) {
		return position % length_ + 0.5;
	}

	double CalculateInitialY(int position) {
		int division = 0;
		if (dimensions_ == 2) {
			division = position / length_;
		}
		else if (dimensions_ == 3) {
			division = (position / length_) % length_;
		}
		return division + 0.5;
	}

	double CalculateInitialZ(int position) {
		int division = 0;
		if (dimensions_ == 3) {
			division = position / pow(length_, 2);
		}
		return division + 0.5;
	}

	double GetX(int particle) {
		return positions_[0][particle];
	}
	double GetY(int particle) {
		return positions_[1][particle];
	}
	double GetZ(int particle) {
		return positions_[2][particle];
	}

public:
	Box(int particles, int length, int dimensions) {
		particles_ = particles;
		length_ = length;
		dimensions_ = dimensions;
		InitialPositions();
	}
	~Box() {
		positions_.~vector();
	}

	void WritePositions(int iteration, double mason, double amplitude_relationship, int repetition, std::string tag) {
		std::ofstream file{ "positions/positions-" + std::to_string(mason) + "-" + std::to_string(amplitude_relationship) + "-" + std::to_string(repetition) + "-" + std::to_string(iteration) + "-" + tag + ".csv" };

		file << "x,y,z\n";

		for (int i = 0; i < particles_; i++) {
			file << GetX(i) << "," << GetY(i) << "," << GetZ(i) << "\n";
		}
		file.close();
	}

	std::vector<double> ReturnX() {
		return positions_[0];
	}

	std::vector<double> ReturnY() {
		return positions_[1];
	}

	std::vector<double> ReturnZ() {
		return positions_[2];
	}

	void SetX(double* x) {
		for (int i = 0; i < particles_; i++) {
			positions_[0][i] = x[i];
		}
	}

	void SetY(double* y) {
		for (int i = 0; i < particles_; i++) {
			positions_[1][i] = y[i];
		}
	}

	void SetZ(double* z) {
		for (int i = 0; i < particles_; i++) {
			positions_[2][i] = z[i];
		}
	}
};