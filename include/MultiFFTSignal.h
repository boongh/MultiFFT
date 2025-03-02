#pragma once
#include <vector>
#include <complex>
#include <type_traits>
#include <iostream>

template <typename T>
struct Signal {
	static_assert(std::is_arithmetic<T>::value || std::is_same<T, std::complex<float>>::value || std::is_same<T, std::complex<double>>::value,
		"T must be an arithmetic type or std::complex<float> or std::complex<double>");

public:

	constexpr static double PI = 3.141592653589793238460;
	constexpr static uint16_t MAX_AMPLITUDE = 32767;

	~Signal() {
		data.clear();
	}


	Signal() {
		this->sampleRate = 0;
		this->data = std::vector<T>();
	}

	Signal(const std::vector<T>& data, int sampleRate, int channel) {
		this->data = data;
		this->sampleRate = sampleRate;
		this->channel = channel;
	}
	Signal(int size, int sampleRate, int channel) {
		this->sampleRate = sampleRate;
		this->data = std::vector<T>(size * channel);
		this->channel = channel;
	}


	int sampleRate;
	int channel;
	std::vector<T>data;



	void PrintSignal() {
		for (int i = 0; i < data.size(); ++i) {
			std::cout << data[i] << "\n";
		}
	}

	/// <summary>
	/// Combine both signals into different channel
	/// File must have same sample rate
	/// </summary>
	/// <param name="signal"></param>
	/// <param name="targetChannel"></param>
	/// <returns></returns>
	static Signal<T>* UpChannel(const Signal& signal1, const Signal& signal2) {
		if (signal1.sampleRate != signal2.sampleRate) {
			std::cerr << "Sample rate is not coherent, unable to combine channels" << "\n";
			throw std::invalid_argument("Signal Incoherent in Signal Upchannel");
		}
		int maxLength = std::max(signal1.data.size(), signal2.data.size());
		int numChannelsOut = signal1.channel + signal2.channel;

		Signal<T> result = Signal<T>(maxLength * numChannelsOut, signal1.sampleRate, numChannelsOut);
	}


	//Combining signal a and b with b shifted forward in time by timeShift seconds
	//If samepleRate is not the same, the signal with the higher sample rate will be used
	static Signal<T> CombineSignal(const Signal& base, const Signal& b, double timeShift, double factor) {

		int sampleratemax = std::max(base.sampleRate, b.sampleRate);
		int maxChannel = std::max(base.channel, b.channel);




		int timeshiftSR = sampleratemax * timeShift;

		int baseSize = base.data.size() * sampleratemax / base.sampleRate;
		int bSize = b.data.size() * sampleratemax / b.sampleRate + sampleratemax * timeShift;

		Signal<T> result = Signal<T>(std::max(baseSize, bSize), sampleratemax, maxChannel);

		

		// Start with adding in base signal at the start
		// Adjust the interval of sample in data to be the same as the sample rate

		// Increment adjust for sample rate difference
		int shift = (static_cast<double>(result.sampleRate) / base.sampleRate);

		for (int i = 0; i < base.data.size(); i++) {
			result.data[i * shift] += base.data[i] * (1 - factor);
		}

		// Calculate shift for signal b;
		shift = (static_cast<double>(result.sampleRate) / b.sampleRate);

		for (int i = 0; i < b.data.size(); i++) {
			result.data[i * shift + timeshiftSR] += b.data[i] * factor;
		}

		return result;
	}

	static Signal<T> GenerateSineWave(int size, int frequency, int sampleRate, double amplitudeNormalized = 1) {

		std::vector<double> sinLookupTable;
		const double phaseIncrement = 2.0 * PI * frequency / sampleRate;

		// Precompute one period of the sine wave
		int samplesPerPeriod = static_cast<int>(std::round(static_cast<double>(sampleRate) / frequency));
		sinLookupTable.reserve(samplesPerPeriod);
		for (int i = 0; i < samplesPerPeriod; ++i) {
			sinLookupTable[i] = (std::sin(i * phaseIncrement));
		}

		Signal<double> result = Signal<double>(size, sampleRate, 1);
		for (int i = 0; i < size; i++) {
			result.data[i] = amplitudeNormalized * sin(2 * PI * frequency * i / (double)sampleRate);
		}
		return (result.data.empty()) ? Signal<T>() : result;
	}

	static Signal<T> ConvertPCMToSignal(Signal<uint16_t> signalIn) {
		Signal<T> dataOut = Signal<T>(static_cast<int>(signalIn.data.size()), signalIn.sampleRate, signalIn.channel);

		for (int i = 0; i < signalIn.data.size(); i++) {
			dataOut.data[i] = static_cast<T>(signalIn.data[i]);
		}
		return (dataOut.data.empty()) ? Signal<T>() : dataOut;
	}

	static Signal<uint16_t> ConvertSignalToPCM(Signal<T> signalIn) {

		Signal<uint16_t> dataOut = Signal<uint16_t>(static_cast<int>(signalIn.data.size()), signalIn.sampleRate, signalIn.channel);

		for (int i = 0; i < signalIn.data.size(); i++) {
			dataOut.data[i] = static_cast<uint16_t>((signalIn.data[i] + 1) * UINT16_MAX);
		}
		return (dataOut.data.empty()) ? Signal<uint16_t>() : dataOut;
	}
};
