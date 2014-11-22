#include "main.h";

using namespace std;

// ���������, ����������� ��������� WAV �����.
struct WAVHEADER
{
	// WAV-������ ���������� � RIFF-���������:

	// �������� ������� "RIFF" � ASCII ���������
	// (0x52494646 � big-endian �������������)
	char chunkId[4];

	// 36 + subchunk2Size, ��� ����� �����:
	// 4 + (8 + subchunk1Size) + (8 + subchunk2Size)
	// ��� ���������� ������ �������, ������� � ���� �������.
	// ����� ������, ��� ������ ����� - 8, �� ����,
	// ��������� ���� chunkId � chunkSize.
	unsigned long chunkSize;

	// �������� ������� "WAVE"
	// (0x57415645 � big-endian �������������)
	char format[4];

	// ������ "WAVE" ������� �� ���� ����������: "fmt " � "data":
	// ���������� "fmt " ��������� ������ �������� ������:

	// �������� ������� "fmt "
	// (0x666d7420 � big-endian �������������)
	char subchunk1Id[4];

	// 16 ��� ������� PCM.
	// ��� ���������� ������ ����������, ������� � ���� �������.
	unsigned long subchunk1Size;

	// ����� ������, ������ ������ ����� �������� ����� http://audiocoding.ru/wav_formats.txt
	// ��� PCM = 1 (�� ����, �������� �����������).
	// ��������, ������������ �� 1, ���������� ��������� ������ ������.
	unsigned short audioFormat;

	// ���������� �������. ���� = 1, ������ = 2 � �.�.
	unsigned short numChannels;

	// ������� �������������. 8000 ��, 44100 �� � �.�.
	unsigned long sampleRate;

	// sampleRate * numChannels * bitsPerSample/8
	unsigned long byteRate;

	// numChannels * bitsPerSample/8
	// ���������� ���� ��� ������ ������, ������� ��� ������.
	unsigned short blockAlign;

	// ��� ���������� "��������" ��� �������� ��������. 8 ���, 16 ��� � �.�.
	unsigned short bitsPerSample;

	// ���������� "data" �������� �����-������ � �� ������.

	// �������� ������� "data"
	// (0x64617461 � big-endian �������������)
	char subchunk2Id[4];

	// numSamples * numChannels * bitsPerSample/8
	// ���������� ���� � ������� ������.
	unsigned long subchunk2Size;

	// ����� ������� ��������������� Wav ������.
};

float bytesToFloat(char firstByte, char secondByte)
{
	unsigned short int s = ((unsigned short int)secondByte << 8) + firstByte;
	float res;
	if (s >= 32768)
		res = (s - 65536) / 32768.0;
	else
		res = s / 32768.0;
	return floor(res * 10000 + 0.5) / 10000;
}

int _tmain(int argc, _TCHAR* argv[])
{
	FILE *file;
	errno_t err;
	err = fopen_s(&file, "obr1.wav", "rb");
	if (err)
	{
		printf_s("Failed open file, error %d", err);
		return 0;
	}

	WAVHEADER header;

	fread_s(&header, sizeof(WAVHEADER), sizeof(WAVHEADER), 1, file);
	vector<float> speech(header.subchunk2Size / 2);
	for (std::vector<float>::iterator it = speech.begin(); it != speech.end(); ++it)
	{
		char a, b;
		fread_s(&a, sizeof(a), sizeof(a), 1, file);
		fread_s(&b, sizeof(b), sizeof(b), 1, file);
		*it = bytesToFloat(a, b);
	}
	

	// ������� ���������� ������
	printf_s("Sample rate: %d\n", header.sampleRate);
	printf_s("Channels: %d\n", header.numChannels);
	printf_s("Bits per sample: %d\n", header.bitsPerSample);
	printf_s("sizeof(header): %d\n", sizeof(short int));

	// ��������� ������������ ��������������� � ��������
	float fDurationSeconds = 1.f * header.subchunk2Size / (header.bitsPerSample / 8) / header.numChannels / header.sampleRate;
	int iDurationMinutes = (int)floor(fDurationSeconds) / 60;
	fDurationSeconds = fDurationSeconds - (iDurationMinutes * 60);
	printf_s("Duration: %02d:%02.f\n", iDurationMinutes, fDurationSeconds);
	for (std::vector<float>::iterator it = speech.begin(); it != speech.end(); ++it)
	{
		cout << *it << endl;
	}
	fclose(file);

	_getch();
	return 0;
}