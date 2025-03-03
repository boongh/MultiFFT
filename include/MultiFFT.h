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

	bool IsPowerOfTwo();
	FrequencyDomain Naive1DDFT(const Signal<double>& input);
	FrequencyDomain FirstIttFFT(const Signal<double>& sample);
	FrequencyDomain SecondIttFFT(const Signal<double>& sample);
	Signal<double> HannWindowing(const Signal<double>& sample);
	bool IsPowerOfTwo(int x);
	bool IsPowerOfTwo(size_t x);
	bool IsPowerOfTwo(unsigned long x);
	static std::vector<std::complex<double>> FirstITTFFTRC(std::vector<std::complex<double>> sample, int n);
	static void SecondITTFFTRC(const std::vector<std::complex<double>>& sample, std::vector<std::complex<double>>& fbins, int start, int stride, int n);
	FrequencyDomain ThridITTFFT(const Signal<double>& sample);
	void ThirdIttFFTRC(std::vector<std::complex<double>>& fbins, const unsigned int SIZE, unsigned int n);
	std::vector<std::complex<double>> ConvertSignalToComplex(const std::vector<double>& input);
	constexpr unsigned int BitReverse(unsigned int b, unsigned int head);
	void PrintMessage(const std::string&);

}