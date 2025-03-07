#pragma once

#include <vector>
#include <complex>


struct FrequencyDomain {
public:
	FrequencyDomain();

	std::vector<std::complex<double>> fbins;
	
	size_t sampleRate;

	FrequencyDomain(size_t bins) {
		fbins = std::vector<std::complex<double>>(bins, 0.0);
		sampleRate = 0;
	}
};
