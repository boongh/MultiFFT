#include "SignalGenerator.h"
#include <vector>
#include "Signal.h"

using namespace SignalGenerator;

Signal<double> GenerateSineWave(int frequency, double amplitude, int sampleRate, int durationSeconds) {

	//Precompute the sin lookup table for fast sin
	std::vector<double> sinLookupTable = std::vector<double>(frequency);
	for (int i = 0; i < frequency; i++) {
		sinLookupTable[i] = sin(2 * PI * i / frequency);
	}


	Signal<double> result = Signal<double>();
	double period = 1.0 / frequency;
	double omega = 2 * PI / period;
	for (int i = 0; i < sampleRate * durationSeconds; i++) {
		result.data[i] = amplitude * sinLookupTable[i];
	}
	return result;
}

//Signal<uint16_t> ConvertSignalToPCM(std::vector<double> data) {
//	Signal<uint16_t> dataOut = Signal<uint16_t>(data.size(), data.);
//	for (int i = 0; i < data.size(); i++) {
//		dataOut[i] = static_cast<uint16_t>((data[i] + 96) * APLITUDE_CONVERSION_FACTOR);
//	}
//
//	return dataOut;
//}