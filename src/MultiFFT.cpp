#include "MultiFFT.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <complex>
#include <math.h>
#include "FrequencyDomain.h"

namespace MultiFFT {

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

		result.sampleRate = samplerate;

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

		result.sampleRate = samplerate;

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

		result.sampleRate = samplerate;

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


	FrequencyDomain ThridITTFFT(const Signal<double>& sample) {
		int N = sample.data.size();

		if (!IsPowerOfTwo(N)) throw std::invalid_argument("Size is not a power of 2");
	
		int log2n = log2(N);

		double temp1;
		double temp2;
		FrequencyDomain result(N);

		//Reverse bit order of sample into freqdomain
		for (int i = 0; i < N; ++i) {
			int cousin = BitReverse(i, log2n);
			if (i < cousin) {
				result.fbins[i] = sample.data[cousin];
				result.fbins[cousin] = sample.data[i];
			}
		}


		//Generate FT
		for (unsigned int n = 2; n < N; n <<= 1) {
			ThirdIttFFTRC(result.fbins, N, n);
		}

		return result;
	}

	/// <summary>
	/// Generate FFT/ Sub ffts
	/// </summary>
	/// <param name="sample">The bit reversed sample</param>
	/// <param name="fbins"></param>
	/// <param name="n"></param>
	/// <param name="stride"></param>
	static void ThirdIttFFTRC(std::vector<std::complex<double>>& fbins, const unsigned int SIZE, unsigned int n) {
		double inversen = 1.0 / (2.0 * n);
		unsigned int stride = n;
		int N = fbins.size();

		for (int k = 0; k < SIZE; k += 2 * n) {
			for (int j = 0; j < stride; j++) {
				std::complex<double> W = std::polar(1.0, -2.0 * PI * j * inversen);

				unsigned int evenIndex = k + j;
				unsigned int oddIndex = k + j + stride;

				std::complex<double> odd = fbins[oddIndex];
				std::complex<double> even = fbins[evenIndex];


				fbins[k + j] = even + W * odd;
				fbins[k + j + stride] = even - W * odd;
			}
		}
	}


	/// <summary>
	/// Convert a running signal of double into complex
	/// </summary>
	/// <param name="input"></param>
	/// <returns></returns>
	std::vector<std::complex<double>> ConvertSignalToComplex(const std::vector<double>& input) {
		int N = input.size();
		std::vector<std::complex<double>> result(N);

		for (int i = 0; i < N; i++) {
			result[i] = std::complex<double>(input[i], 0.0);
		}

		return result;
	}

	constexpr unsigned int BitReverse(unsigned int b, unsigned int head) {

		b = (b & 0xFFFF0000) >> 16 | (b & 0x0000FFFF) << 16;
		b = (b & 0xFF00FF00) >> 8 | (b & 0x00FF00FF) << 8;
		b = (b & 0xF0F0F0F0) >> 4 | (b & 0x0F0F0F0F) << 4;
		b = (b & 0xCCCCCCCC) >> 2 | (b & 0x33333333) << 2;
		b = (b & 0xAAAAAAAA) >> 1 | (b & 0x55555555) << 1;

		b >>= (sizeof(unsigned int) * 8 - head);
		return b;
	}


	Signal<double> HannWindowing(const Signal<double>& sample) {
		int N = sample.data.size();
		Signal<double> windowedSignal(N, sample.sampleRate, sample.channel);

		for (int i = 0; i < N; ++i) {
			windowedSignal.data[i] = sample.data[i] * sin(PI * i / N);
		}

		return windowedSignal;
	}

	bool IsPowerOfTwo(int x)
	{
		return (x != 0) && ((x & (x - 1)) == 0);
	}

	bool IsPowerOfTwo(size_t x)
	{
		return (x != 0) && ((x & (x - 1)) == 0);
	}

	bool IsPowerOfTwo(unsigned long x)
	{
		return (x != 0) && ((x & (x - 1)) == 0);
	}
}