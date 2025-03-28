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

	/// <summary>
	/// Naive implementation of the DFT
	/// </summary>
	/// <param name="input"></param>
	/// <returns>Frequency domain of the signal</returns>
	FrequencyDomain Naive1DDFT(const Signal<double>& input);

	/// <summary>
	/// First itteration of the FFT
	/// </summary>
	/// <param name="sample"></param>
	/// <returns></returns>
	FrequencyDomain FirstIttFFT(const Signal<double>& sample);
	static std::vector<std::complex<double>> FirstITTFFTRC(std::vector<std::complex<double>> sample, int n);

	FrequencyDomain SecondIttFFT(const Signal<double>& sample);
	static void SecondITTFFTRC(const std::vector<std::complex<double>>& sample, std::vector<std::complex<double>>& fbins, int start, int stride, int n);

	FrequencyDomain ThridITTFFT(const Signal<double>& sample);
	void ThirdIttFFTRC(std::vector<std::complex<double>>& fbins);

	FrequencyDomain FourthIttFTT(const Signal<double>& sample);
	static void FourthITTFTTInternal(std::vector<std::complex<double>>& fbins);

	FrequencyDomain FifthIttFFT(const Signal<double>& sample);
	static void FifthFFTInternal(std::vector<std::complex<double>>& fbins);

	FrequencyDomain SixthIttFFT(const Signal<double>& sample);
	static void SixthFFTInternal(std::vector<std::complex<double>>& fbins);

	/// <summary>
	/// Inverse FFT by inverting the SixthIttFFT
	/// </summary>
	/// <param name="FD"></param>
	/// <returns>New Signal</returns>
	Signal<std::complex<double>> InverseFFT(FrequencyDomain FD);

	std::vector<std::complex<double>> ConvertSignalToComplex(const std::vector<double>& input);

	unsigned int RoundDownPowerOfTwo(unsigned int n);

	inline constexpr unsigned int BitReverse(unsigned int b, unsigned int head);

	Signal<double> HannWindowing(const Signal<double>& sample);

	void PrintMessage(const std::string&);

	bool IsPowerOfTwo(int x);
	bool IsPowerOfTwo(size_t x);
	bool IsPowerOfTwo(unsigned long x);

}