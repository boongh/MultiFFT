#include "Signal.h"



//template<typename T>
//Signal<T>::Signal(const vector<T>& data, int sampleRate, double amplitude){
//	this->data = data;
//	this->sampleRate = sampleRate;
//	this->amplitude = amplitude;
//}
//
//template<typename T>
//Signal<T>::Signal(int size, int sampleRate, double amplitude){
//	this->sampleRate = sampleRate;
//	this->amplitude = amplitude;
//	this->data = vector<T>(size);
//}
//
//
//template<typename T>
//Signal<T> Signal<T>::CombineSignal(Signal<T> base, Signal<T> b, double timeShift){
//	Signal<T> result = Signal<T>(max(base.data.size(), b.data.size()), max(base.sampleRate, b.sampleRate), base.amplitude);
//
//	// Expand the signal range as required if signal b is longer considering the shift and their sample rate
//	result.sampleRate = max(base.sampleRate, b.sampleRate);
//	int baseSize = base.data.size() * result.sampleRate / base.sampleRate;
//	int bSize = b.data.size() * result.sampleRate / b.sampleRate + result.sampleRate * timeShift;
//	result.data.resize(max(baseSize, bSize));
//
//	// Start with adding in base signal at the start
//	// Adjust the interval of sample in data to be the same as the sample rate
//
//	// Increment adjust for sample rate difference
//	int shift = (result.sampleRate / base.sampleRate);
//	int shiftedSample;
//
//	for (int i = 0; i < base.data.size(); i ++) {
//		shiftedSample = i * result.sampleRate / base.sampleRate;
//		result.data[shiftedSample] = base.data[i];
//	}
//
//	// Calculate shift for signal b;
//	shift = (result.sampleRate / b.sampleRate);
//
//	for (int i = 0; i < b.data.size(); i++) {
//		shiftedSample = i * result.sampleRate / b.sampleRate + result.sampleRate * timeShift;
//		result.data[shiftedSample] += b.data[i];
//	}
//
//	return result;
//}
//
//
//template<typename T>
//Signal<T> Signal<T>::GenerateSineWave(int size, int frequency, int sampleRate, int amplitude){
//
//	std::vector<double> sinLookupTable;
//	const double phaseIncrement = 2.0 * PI * frequency / sampleRate;
//
//	// Precompute one period of the sine wave
//	int samplesPerPeriod = static_cast<int>(std::round(static_cast<double>(sampleRate) / frequency));
//	sinLookupTable.reserve(samplesPerPeriod);
//	for (int i = 0; i < samplesPerPeriod; ++i) {
//		sinLookupTable[i] = (std::sin(i* phaseIncrement));
//	}
//
//	Signal<double> result = Signal<double>(size);
//	for (int i = 0; i < size; i++) {
//		result.data[i] = amplitude * sinLookupTable[i % samplesPerPeriod];
//	}
//	return (result.data.empty()) ? Signal<T>() : result;
//}
//
//template<typename T>
//Signal<uint16_t> Signal<T>::ConvertSignalToPCM(Signal<T> signalIn, double amplitude)
//{
//	Signal<uint16_t> dataOut = Signal<uint16_t>(signalIn.data.size(), signalIn.sampleRate, signalIn.amplitude);
//
//	for (int i = 0; i < signalIn.data.size(); i++) {
//		dataOut.data[i] = static_cast<uint16_t>((signalIn.data[i] + amplitude) * APLITUDE_CONVERSION_FACTOR);
//	}
//	return (dataOut.data.empty) ? Signal<uint16_t>() : dataOut;
//}
