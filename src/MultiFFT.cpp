#include "include/MultiFFT.h" as MultiFFT
#include <iostream>
#include <math.h>
#include <vector>
#include <complex>


using namespace std;

namespace MultiFFT {
	void PrintMessage(const std::string& str)
	{
		std::cout << str << std::endl;
	}

	vector<complex<double>> Naive1DDFT(const vector<int>& input) {
		vector<complex<double>> result;
		const double PI = 3.141592653589793238460;
		const complex<double> i(0, 1);
		int N = input.size();
		int halfN = N / 2;


		//Loop through the input vector


		return result;
	}
}