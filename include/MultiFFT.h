#pragma once

#include <string>
#include <vector>
#include <complex>

using namespace std;

namespace MultiFFT {
	void PrintMessage(const std::string&);
	vector<complex<double>> Naive1DDFT(const vector<int>&);
}