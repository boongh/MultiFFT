#pragma once

#include <vector>
#include <complex>


class FrequencyDomain {
public:
	std::vector<std::complex<double>> fbins;
	double binDifFreq;
	int startBin;

	FrequencyDomain(int bins) {
		fbins = std::vector<std::complex<double>>(bins, 0.0);
		binDifFreq = 0.0;
		startBin = 0;
	}
};
