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
		size_t N = input.data.size();
		size_t halfN = N / 2;
		double inverseN = 1.0 / N;
		size_t samplerate = input.sampleRate;

		FrequencyDomain result(N);

		result.sampleRate = samplerate;

		int hundredThousand = 1000;
		int counter = 0;

		for (size_t i = 0; i < N; ++i) {
			std::complex<double> angle = -IUNIT * (2.0 * PI * i * inverseN);
			for (size_t k = 0; k < N; ++k) {
				result.fbins[i] += input.data[k] * exp(angle * (double)k) / static_cast<double>(halfN);
			}
		}


		//Loop through the input vector


		return result;
	}


	/// <summary>
	/// First itteration of the FFT implementation
	/// Sample has to be power of 2
	/// Very big memory overhead
	/// </summary>
	/// <param name="sample"></param>
	/// <returns></returns>
	FrequencyDomain FirstIttFFT(const Signal<double>& sample) {
		size_t N = static_cast<unsigned long long>(sample.data.size());

		if (!IsPowerOfTwo(N) && N != 1) {
			throw std::invalid_argument("n not power of 2 FFT");
		}

		size_t halfN = N / 2;
		size_t samplerate = sample.sampleRate;

		FrequencyDomain result(N);

		result.sampleRate = samplerate;

		std::vector<std::complex<double>> conversionToComplex(N);

		for (size_t i = 0; i < N; i++) {
			conversionToComplex[i] = std::complex<double>(sample.data[i], 0.0);
		};

		result.fbins = FirstITTFFTRC(conversionToComplex, N);

		for (int i = 0; i < N; i++) {
			result.fbins[i] /= static_cast<double>(halfN);
		}

		return result;
	}


	/// <summary>
	/// Internal helper function for the first itteration of FFT
	/// </summary>
	/// <param name="sample"></param>
	/// <param name="n"></param>
	/// <returns></returns>
	static std::vector<std::complex<double>> FirstITTFFTRC(std::vector<std::complex<double>> sample, int n) {

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

		if (!IsPowerOfTwo(N) && N != 1) {
			throw std::invalid_argument("n not power of 2 FFT");
		}

		int halfN = N / 2;
		int samplerate = sample.sampleRate;

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



	/// <summary>
	/// Interface for FFT
	/// Third itterationf FFT algorithm
	/// Bit reversal permutation + one allocation
	/// No recursion, NASA likes :)
	/// </summary>
	/// <param name="sample"></param>
	/// <returns></returns>
	/// 
	FrequencyDomain ThridITTFFT(const Signal<double>& sample) {
		size_t N = sample.data.size();

		if (!IsPowerOfTwo(N)) throw std::invalid_argument("Size is not a power of 2");
	
		size_t log2n = static_cast<size_t>(log2(N));

		FrequencyDomain result(N);

		//Reverse bit order of sample into freqdomain
		for (size_t i = 0; i < N; ++i) {
			size_t cousin = BitReverse(i, log2n);
			result.fbins[i] = sample.data[cousin];
			result.fbins[cousin] = sample.data[i];
		}
		ThirdIttFFTRC(result.fbins);

		return result;
	}

	/// <summary>
	/// Generate FFT/ Sub ffts
	/// </summary>
	/// <param name="sample">The bit reversed sample</param>
	/// <param name="fbins"></param>
	/// <param name="n"></param>
	/// <param name="stride"></param>
	static void ThirdIttFFTRC(std::vector<std::complex<double>>& fbins) {
		size_t SIZE = fbins.size();

		for (size_t n = 2; n < SIZE; n <<= 1) {

			std::cout << "First layer " << n << "\n";
			double inversen = 0.5 / n;
			size_t stride = n;

			for (size_t k = 0; k < SIZE; k += 2 * n) { //Completes every butterfly loop

				for (size_t j = 0; j < stride; j++) { //Completes each butterfly

					size_t evenIndex = k + j;
					size_t oddIndex = k + j + stride;

					std::complex<double> odd = fbins[oddIndex];
					std::complex<double> even = fbins[evenIndex];


					//std::complex<double> W = std::polar(1.0, -2.0 * PI * j * inversen);
					//std::complex<double> odd = fbins[oddIndex];

					std::complex<double> W = std::polar(1.0, -2.0 * PI * j * inversen) * odd;

					fbins[k + j] = even + W;
					fbins[k + j + stride] = even - W;
				}
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

	unsigned int RoundDownPowerOfTwo(unsigned int n) {
		return static_cast<int>(std::log2(n));
	}

	inline constexpr unsigned int BitReverse(unsigned int b, unsigned int head) {

		b = (b & 0xFFFF0000) >> 16 | (b & 0x0000FFFF) << 16;
		b = (b & 0xFF00FF00) >> 8 | (b & 0x00FF00FF) << 8;
		b = (b & 0xF0F0F0F0) >> 4 | (b & 0x0F0F0F0F) << 4;
		b = (b & 0xCCCCCCCC) >> 2 | (b & 0x33333333) << 2;
		b = (b & 0xAAAAAAAA) >> 1 | (b & 0x55555555) << 1;

		b >>= (sizeof(unsigned int) * 8 - head);
		return b;
	}


	Signal<double> HannWindowing(const Signal<double>& sample) {
		unsigned long long N = static_cast<unsigned long long>(sample.data.size());
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