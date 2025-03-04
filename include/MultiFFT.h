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
	FrequencyDomain ThridITTFFT(const Signal<double>& sample);
	static std::vector<std::complex<double>> FirstITTFFTRC(std::vector<std::complex<double>> sample, int n);
	static void SecondITTFFTRC(const std::vector<std::complex<double>>& sample, std::vector<std::complex<double>>& fbins, int start, int stride, int n);
	void ThirdIttFFTRC(std::vector<std::complex<double>>& fbins);
	std::vector<std::complex<double>> ConvertSignalToComplex(const std::vector<double>& input);
	unsigned int RoundDownPowerOfTwo(unsigned int n);
	constexpr unsigned long long BitReverse(unsigned long long b, unsigned long long head);
	Signal<double> HannWindowing(const Signal<double>& sample);
	void PrintMessage(const std::string&);
	bool IsPowerOfTwo(int x);
	bool IsPowerOfTwo(size_t x);
	bool IsPowerOfTwo(unsigned long x);

}