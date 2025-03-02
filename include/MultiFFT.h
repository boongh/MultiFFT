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
	FrequencyDomain SecondIttFFT(const Signal<double>& sample);
	Signal<double> HannWindowing(const Signal<double>& sample);
	static std::vector<std::complex<double>> FirstITTFFTRC(std::vector<std::complex<double>> sample, int n);
	static void SecondITTFFTRC(const std::vector<std::complex<double>>& sample, std::vector<std::complex<double>>& fbins, int start, int stride, int n);
	int BitReverse(int x, int log2n);
	void PrintMessage(const std::string&);

}