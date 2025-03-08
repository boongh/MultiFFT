#pragma once

#include <vector>
#include <complex>


struct FrequencyDomain {
public:
	FrequencyDomain();

	alignas(32) std::vector<std::complex<double>> fbins;
	
	size_t sampleRate;

	FrequencyDomain(size_t bins) {
		fbins = std::vector<std::complex<double>>(bins);
		sampleRate = 0;
	}
};
