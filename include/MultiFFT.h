#pragma once

#include <string>
#include <vector>
#include <complex>
#include "MultiFFTSignal.h"
#include "FrequencyDomain.h"

using namespace std;

namespace MultiFFT {
	constexpr std::complex<double> IUNIT(0.0, 1.0);
	constexpr double PI = 3.141592653589793238460;

	FrequencyDomain Naive1DDFT(const Signal<double>& input);
	FrequencyDomain FirstIttFFT(const Signal<double>& sample);
	static std::vector<std::complex<double>> FirstITTFFTRC(std::vector<std::complex<double>> sample, int start, int stride, int n);
	void PrintMessage(const std::string&);

}