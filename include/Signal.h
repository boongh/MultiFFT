#pragma once
#include <vector>
#include <complex>
#include <type_traits>

using namespace std;


template <typename T>
class Signal {
	static_assert(is_arithmetic<T>::value || is_same<T, complex<float>>::value || is_same<T, complex<double>>::value,
		"T must be an arithmetic type or std::complex<float> or std::complex<double>");

public:

	constexpr static double PI = 3.141592653589793238460;
	constexpr static uint16_t MAX_AMPLITUDE = 32767;
	constexpr static uint16_t MAX_AMPLITUDE_DB = 96;
	constexpr static uint16_t AMPLITUDE_CONVERSION_FACTOR = MAX_AMPLITUDE / MAX_AMPLITUDE_DB;

	~Signal() {
		data.clear();
	}


	Signal() {
		this->sampleRate = 0;
		this->amplitude = 0;
		this->data = vector<T>();
	}

	Signal(const vector<T>& data, int sampleRate, double amplitude, int channel) {
		this->data = data;
		this->sampleRate = sampleRate;
		this->amplitude = amplitude;
		this->channel = channel;
	}
	Signal(int size, int sampleRate, double amplitude, int channel) {
		this->sampleRate = sampleRate;
		this->amplitude = amplitude;
		this->data = vector<T>(size);
		this->channel = channel;
	}


	int sampleRate;
	double amplitude;
	int channel;
	vector<T> data;



	//Combining signal a and b with b shifted forward in time by timeShift seconds
	//If samepleRate is not the same, the signal with the higher sample rate will be used
	static Signal<T> CombineSignal(Signal base, Signal b, double timeShift, double factor) {

		int sampleratemax = max(base.sampleRate, b.sampleRate);

		int maxChannel = max(base.channel, b.channel);

		int timeshiftSR = sampleratemax * timeShift;

		int baseSize = base.data.size() * sampleratemax / base.sampleRate;
		int bSize = b.data.size() * sampleratemax / b.sampleRate + sampleratemax * timeShift;

		Signal<T> result = Signal<T>(max(baseSize, bSize), sampleratemax, base.amplitude, maxChannel);

		

		// Start with adding in base signal at the start
		// Adjust the interval of sample in data to be the same as the sample rate

		// Increment adjust for sample rate difference
		int shift = (static_cast<double>(result.sampleRate) / base.sampleRate);

		for (int i = 0; i < base.data.size(); i++) {
			result.data[i * shift] += base.data[i];
		}

		// Calculate shift for signal b;
		shift = (static_cast<double>(result.sampleRate) / b.sampleRate);

		for (int i = 0; i < b.data.size(); i++) {
			result.data[i * shift + timeshiftSR] += b.data[i];
		}

		return result;
	}

	static Signal<T> GenerateSineWave(int size, int frequency, int sampleRate, int amplitudeDB = MAX_AMPLITUDE_DB) {

		std::vector<double> sinLookupTable;
		const double phaseIncrement = 2.0 * PI * frequency / sampleRate;

		// Precompute one period of the sine wave
		int samplesPerPeriod = static_cast<int>(std::round(static_cast<double>(sampleRate) / frequency));
		sinLookupTable.reserve(samplesPerPeriod);
		for (int i = 0; i < samplesPerPeriod; ++i) {
			sinLookupTable[i] = (std::sin(i * phaseIncrement));
		}

		Signal<double> result = Signal<double>(size, sampleRate, amplitudeDB, 1);
		for (int i = 0; i < size; i++) {
			result.data[i] = amplitudeDB * sinLookupTable[i % samplesPerPeriod];
		}
		return (result.data.empty()) ? Signal<T>() : result;
	}

	static Signal<T> ConvertPCMToSignal(Signal<uint16_t> signalIn, double amplitude = MAX_AMPLITUDE_DB) {
		Signal<T> dataOut = Signal<T>(static_cast<int>(signalIn.data.size()), signalIn.sampleRate, signalIn.amplitude, signalIn.channel);

		for (int i = 0; i < signalIn.data.size(); i++) {
			dataOut.data[i] = static_cast<T>(signalIn.data[i]) * amplitude / 32768.0;
		}
		return (dataOut.data.empty()) ? Signal<T>() : dataOut;
	}

	static Signal<uint16_t> ConvertSignalToPCM(Signal<T> signalIn, double amplitude = MAX_AMPLITUDE_DB) {
		double conversionFactor = AMPLITUDE_CONVERSION_FACTOR * amplitude / MAX_AMPLITUDE;

		Signal<uint16_t> dataOut = Signal<uint16_t>(static_cast<int>(signalIn.data.size()), signalIn.sampleRate, signalIn.amplitude, signalIn.channel);

		for (int i = 0; i < signalIn.data.size(); i++) {
			dataOut.data[i] = static_cast<uint16_t>((signalIn.data[i] + amplitude) * AMPLITUDE_CONVERSION_FACTOR);
		}
		return (dataOut.data.empty()) ? Signal<uint16_t>() : dataOut;
	}
};
