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

	Signal(const std::vector<T>& data, size_t sampleRate, size_t channel) {
		this->data = data;
		this->sampleRate = sampleRate;
		this->channel = channel;
	}
	Signal(size_t size, size_t sampleRate, size_t channel) {
		this->sampleRate = sampleRate;
		this->data = std::vector<T>(size * channel);
		this->channel = channel;
	}


	size_t sampleRate;
	size_t channel;
	std::vector<T>data;



	void PrintSignal() {
		for (size_t i = 0; i < data.size(); ++i) {
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
	static Signal<T> UpChannelCombine(const Signal& signal1, const Signal& signal2) {
		if (signal1.sampleRate != signal2.sampleRate) {
			std::cerr << "Sample rate is not coherent, unable to combine channels" << "\n";
			throw std::invalid_argument("Sample rate  Incoherent in Signal Upchannel");
		}
		// Calculate the number of samples per channel for each signal
		int samples1 = signal1.data.size() / signal1.channel;
		int samples2 = signal2.data.size() / signal2.channel;

		int maxLength = std::max(samples1, samples2);
		int numChannelsOut = signal1.channel + signal2.channel;

		Signal<T> result = Signal<T>(maxLength, signal1.sampleRate, numChannelsOut);

		// For each sample, copy available channels from signal1 and signal2 into the result.
		for (int i = 0; i < maxLength; i++) {
			int resultBase = i * numChannelsOut;
			// For signal1 channels
			for (int j = 0; j < signal1.channel; j++) {
				// If signal1 has this sample, copy it; otherwise, fill with zero.
				if (i < samples1)
					result.data[resultBase + j] = signal1.data[i * signal1.channel + j];
				else
					result.data[resultBase + j] = static_cast<T>(0);
			}

			// For signal2 channels
			for (int j = 0; j < signal2.channel; j++) {
				if (i < samples2)
					//Offset after signal 1
					result.data[resultBase + signal1.channel + j] = signal2.data[i * signal2.channel + j];
				else
					result.data[resultBase + signal1.channel + j] = static_cast<T>(0);
			}
		}

		return result;
	}


	/// <summary>
	/// Split up all the channels a single signal into their own signals
	/// </summary>
	/// <param name="signalin"></param>
	/// <returns>vector of signals</returns>
    static std::vector<Signal<T>> SplitChannels(const Signal& signalin) {
		int channels = signalin.channel;
		int newSize = signalin.data.size() / channels;

		std::vector<Signal<T>> results(channels, Signal<T>(newSize, signalin.sampleRate, 1));

		for (int c = 0; c < channels; c++) {
			for (int j = 0; j < newSize; ++j) {
				results[c].data[j] = signalin.data[j * channels + c];
			}
		}

		std::cout << "Split channel" << "\n";
		std::cout << "Signal 1 : " << signalin.data.size() << "\n";
		std::cout << "Datapoint : " << results[0].data.size() << "\n";
		std::cout << "Size of a point : " << sizeof(results[0].data[0]) << "\n";
		std::cout << "Sample rate : " << results[0].sampleRate << "\n";
		std::cout << "Channels : " << results[0].channel << "\n";


		return results;
    }



	//Combining signal a and b with b shifted forward in time by timeShift seconds
	//If samepleRate is not the same, the signal with the higher sample rate will be used
	static Signal<T> CombineSignal(const Signal& base, const Signal& b, double timeShift, double factor) {
		if (base.sampleRate != b.sampleRate) {
			throw std::invalid_argument("Sample rate or channel not coherent");
		}

		int TimeShiftSamples = base.sampleRate * timeShift;

		size_t baseSize = base.data.size() / base.channel;
		size_t bSize = b.data.size() / b.channel;

		size_t maxSize = std::max(baseSize, bSize + TimeShiftSamples);
		int maxChannel = std::max(b.channel, base.channel);

		Signal<T> result = Signal<T>(maxSize, base.sampleRate, maxChannel);
		

		//Base signal
		for (size_t i = 0; i < baseSize; i++) {

			//Insert each channel per smaple
			for (int c = 0; c < base.channel; c++) {
				result.data[i * maxChannel + c] += base.data[i * base.channel + c] * (1 - factor);
			}
		}


		//B signal
		for (size_t i = 0; i < bSize; i++) {

			//Same but with shift
			for (int c = 0; c < b.channel; c++) {
				result.data[(i + TimeShiftSamples) * maxChannel + c] += b.data[i * b.channel + c] * factor;
			}
		}

		/*std::cout << "Combined" << "\n";
		std::cout << "Datapoint : " << result.data.size() << "\n";
		std::cout << "Size of a point : " << sizeof(result.data[0]) << "\n";
		std::cout << "Sample rate : " << result.sampleRate << "\n";
		std::cout << "Channels : " << result.channel<< "\n";*/

		return result;
	}

	static Signal<T> GenerateSineWave(size_t size, int frequency, int sampleRate, double amplitudeNormalized = 1) {

		std::vector<double> sinLookupTable;
		const double phaseIncrement = 2.0 * PI * frequency / sampleRate;

		// Precompute one period of the sine wave
		int samplesPerPeriod = static_cast<int>(std::round(static_cast<double>(sampleRate) / frequency));
		sinLookupTable.resize(samplesPerPeriod);

		Signal<double> result = Signal<double>(size, sampleRate, 1);
		for (size_t i = 0; i < size; i++) {
			result.data[i] = amplitudeNormalized * sin(2 * PI * frequency * i / (double)sampleRate);
		}
		return (result.data.empty()) ? Signal<T>() : result;
	}


	/// <summary>
	/// Convert Signed PCM Signal to normalized signal format
	/// </summary>
	/// <param name="signalIn"></param>
	/// <returns>New normalized signal</returns>
	static Signal<T> ConvertSignedPCMToSignal(const Signal<int16_t>& signalIn) {
		size_t sample = signalIn.data.size() / signalIn.channel;
		Signal<T> dataOut = Signal<T>(sample, signalIn.sampleRate, signalIn.channel);

		for (int i = 0; i < sample; i++) {
			dataOut.data[i] = static_cast<T>(signalIn.data[i]) / 32768; // Signal deals with assumption that signal is normalised and signed
		}
		return (dataOut.data.empty()) ? Signal<T>() : dataOut;
	}


	/// <summary>
	/// Convert signed normalized Signal format into signed PCM
	/// </summary>
	/// <param name="signalIn"></param>
	/// <returns>New scaled PCM signal</returns>
	static Signal<int16_t> ConvertSignalToPCM(const Signal<T>& signalIn) {

		Signal<int16_t> dataOut = Signal<int16_t>(signalIn.data.size() / signalIn.channel, signalIn.sampleRate, signalIn.channel);

		for (int i = 0; i < signalIn.data.size(); i++) {
			dataOut.data[i] = static_cast<uint16_t>(signalIn.data[i] * (-INT16_MIN));
		}

		return (dataOut.data.empty()) ? Signal<int16_t>() : dataOut;
	}

};
