#pragma once

struct WAVHeader {
	char chunkID[4];
	unsigned int chunkSize;
	char format[4];
	char subchunk1ID[4];
	unsigned int subchunk1Size;
	unsigned short audioFormat;
	unsigned short numChannels;
	unsigned int sampleRate;
	unsigned int byteRate;
	unsigned short blockAlign;
	unsigned short bitsPerSample;
	char subchunk2ID[4];
	unsigned int subchunk2Size;
};

struct