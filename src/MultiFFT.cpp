#include "MultiFFT.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <complex>
#include <math.h>
#include "FrequencyDomain.h"

namespace MultiFFT {

	bool IsPowerOfTwo(int x)
	{
		return (x != 0) && ((x & (x - 1)) == 0);
	}
	bool IsPowerOfTwo(unsigned long x)
	{
		return (x != 0) && ((x & (x - 1)) == 0);
	}

	void PrintMessage(const std::string& str)
	{
		std::cout << str << std::endl;
	}

	FrequencyDomain Naive1DDFT(const Signal<double>& input) {
		int N = input.data.size();
		int halfN = N / 2;
		double inverseN = 1.0 / N;
		int samplerate = input.sampleRate;
		int maxF = samplerate / 2;
		int minBin = N / static_cast<double>(samplerate);
		double frequencyResolution = 1.0 / minBin;
		int maxBin = maxF / frequencyResolution;

		FrequencyDomain result(N);

		result.binDifFreq = frequencyResolution;
		result.startBin = 0;

		int hundredThousand = 1000;
		int counter = 0;

		for (int i = 0; i < N; ++i) {
			std::complex<double> angle = -IUNIT * (2.0 * PI * i * inverseN);
			for (int k = 0; k < N; ++k) {
				result.fbins[i] += input.data[k] * exp(angle * (double)k) / static_cast<double>(halfN);
			}

			if (counter++ == hundredThousand) {
				std::cout << result.fbins[i] << "\n";
				std::cout << i << " bin" << "\n";
				counter = 0;
			}
		}


		//Loop through the input vector


		return result;
	}

	FrequencyDomain FirstIttFFT(const Signal<double>& sample) {
		int N = sample.data.size();
		int halfN = N / 2;
		double inverseN = 1.0 / N;
		int samplerate = sample.sampleRate;
		int maxF = samplerate / 2;
		int minBin = N / static_cast<double>(samplerate);
		double frequencyResolution = 1.0 / minBin;
		int maxBin = maxF / frequencyResolution;

		FrequencyDomain result(N);

		result.binDifFreq = frequencyResolution;
		result.startBin = 0;

		std::vector<std::complex<double>> conversionToComplex(N);

		for (int i = 0; i < N; i++) {
			conversionToComplex[i] = std::complex<double>(sample.data[i], 0.0);
		};

		result.fbins = FirstITTFFTRC(conversionToComplex, N);

		for (int i = 0; i < N; i++) {
			result.fbins[i] /= static_cast<double>(halfN);
		}

		return result;
	}

	static std::vector<std::complex<double>> FirstITTFFTRC(std::vector<std::complex<double>> sample, int n) {

		if (!IsPowerOfTwo(n) && n != 1){
			std::invalid_argument("n not power of 2 FFT");
		}

		if (n == 1) { return sample; }
			
		int M = n / 2;
		std::vector<std::complex<double>> freqBin(n, 0);

		std::vector<std::complex<double>> Xeven(M, 0);
		std::vector<std::complex<double>> Xodd(M, 0);

		for (int i = 0; i < M; i++) {
			Xeven[i] = sample[2 * i];
			Xodd[i] = sample[2 * i + 1];
		}

		std::vector<std::complex<double>> Feven(M, 0);
		Feven = FirstITTFFTRC(Xeven, M);
		std::vector<std::complex<double>> Fodd(M, 0);
		Fodd = FirstITTFFTRC(Xodd, M);

		for (int k = 0; k < M; k++) {	
			std::complex<double> W = std::polar(1.0, -2.0 * PI * k / n) * Fodd[k];
			freqBin[k] = (Feven[k] + W);
			freqBin[k + M]	= Feven[k] - W;
		}
		return freqBin;
	}

	FrequencyDomain SecondIttFFT(const Signal<double>& sample) {
		int N = sample.data.size();
		int halfN = N / 2;
		double inverseN = 1.0 / N;
		int samplerate = sample.sampleRate;
		int maxF = samplerate / 2;
		int minBin = N / static_cast<double>(samplerate);
		double frequencyResolution = 1.0 / minBin;
		int maxBin = maxF / frequencyResolution;

		FrequencyDomain result(N);

		result.binDifFreq = frequencyResolution;
		result.startBin = 0;

		std::vector<std::complex<double>> conversionToComplex(N);

		for (int i = 0; i != N - 1; i++) {
			conversionToComplex[i] = std::complex<double>(sample.data[i], 0.0);
		}

		SecondITTFFTRC(conversionToComplex, result.fbins, 0, 1, N);

		for (int i = 0; i < N; i++) {
			result.fbins[i] /= static_cast<double>(halfN);
		}

		return result;
	}

	Signal<double> HannWindowing(const Signal<double>& sample) {
		int N = sample.data.size();
		Signal<double> windowedSignal(N, sample.sampleRate, sample.channel);

		for (int i = 0; i < N; ++i) {
			windowedSignal.data[i] = sample.data[i] * sin(PI * i / N);
		}

		return windowedSignal;
	}

	static void SecondITTFFTRC(const std::vector<std::complex<double>>& sample, std::vector<std::complex<double>>& fbins,  int start, int stride, int n) {

		if (!IsPowerOfTwo(n) && n != 1) {
			std::invalid_argument("n not power of 2 FFT");
		}

		if (n == 1) {
			fbins[0] = sample[start];
			return;
		}

		int M = n / 2;

		std::vector<std::complex<double>> fbinsEven(M);
		std::vector<std::complex<double>> fbinsOdd(M);

		SecondITTFFTRC(sample, fbinsEven, start, stride * 2, M);			//Event indices
		SecondITTFFTRC(sample, fbinsOdd, start + stride, stride * 2, M);	//Odd indices

		for (int k = 0; k < M; k++) {
			std::complex<double> W = std::polar(1.0, -2.0 * PI * k / n);
			fbins[k] = fbinsEven[k] + W * fbinsOdd[k];
			fbins[k + M] = fbinsEven[k] - W * fbinsOdd[k];
		}
	}	

	int BitReverse(int x, int log2n) {
		int reversed = 0;
		for (int i = 0; i < log2n; i++) {
			reversed = (reversed << 1) | (x & 1);
			x >>= 1;
		}
		return reversed;
	}
}