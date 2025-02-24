#pragma once

#include <string>

namespace MultiFFT {
	class MultiFFT {
	public:
		MultiFFT();
		~MultiFFT();
		void fft();
	};

	void PrintMessage(const std::string&);
}