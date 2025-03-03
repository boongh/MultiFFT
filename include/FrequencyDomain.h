#pragma once

#include <vector>
#include <complex>


class FrequencyDomain {
public:
	FrequencyDomain();

	std::vector<std::complex<double>> fbins;
	
	int sampleRate;

	FrequencyDomain(int bins) {
		fbins = std::vector<std::complex<double>>(bins, 0.0);
		sampleRate = 0;
	}
};
