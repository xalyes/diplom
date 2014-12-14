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

void filter(const float b, const float a, const vector<float> &x, vector<float> &speech)
{
	speech[0] = b*x[0];
	float sum1, sum2;
	int N = x.size();
	for (int i = 1; i < N; ++i)
	{
		sum1 = 0;
		sum2 = 0;
		for (int j = i - i*b; j <= i; ++j)
		{
			sum1 += x[j];
		}
		sum1 *= b;
		for (int j = i - i*a; j <= i - 1; ++j)
		{
			sum2 += speech[j];
		}
		sum2 *= a;
		speech[i] = sum1 - sum2;
	}
}

void vec2frames(vector<float> &vec, int Nw, int Ns, vector <vector <float> > &frames)
{
	int L = vec.size();

	int M = floor((L - Nw) / Ns + 1);
	int E = (L - ((M - 1) * Ns + Nw));
	if (E > 0)
	{
		int P = Nw - E;
		vec.insert(vec.end(), P, 0);
		M = M + 1;
	}

	vector<int> indf(M);
	vector<int> inds(Nw);
	for (int i = 0; i < M; ++i)
		indf[i] = i * Ns;
	for (int i = 0; i < Nw; ++i)
		inds[i] = i + 1;

	//vector <vector <float> > frames(inds.size(), vector<float>(indf.size()));
	for (int i = 0; i < inds.size(); ++i)
		for (int j = 0; j < indf.size(); ++j)
		{
			frames[i][j] = vec[inds[i] + indf[j]];
		}


}

bool get_feature_vector(vector<float> speech, const int fs, vector <vector <float>> &MFCCs_train)
{
	bool use_ceplifter = false;
	float alpha = 0.97;				// ����������� ��������������� ���������
	int frame_duration = 100;		// ������������ ������ � ��
	int frame_shift = 50;			// ����� ������� � ��
	int M = 40;						// ���������� ��������
	int L = 22;						// �������� ������������
	float ERBScaleFactor = 0.75;	// ����������� ������ ��������
	int fmin = 100;					// ����������� �������, ��
	int fmax = 6250;				// ������������ �������, ��
	int Ncc = 29;					// ���������� ��� - ������������ �������������

	//filter(1 - alpha, 1, speech, speech);

	int Nw = round(0.001 * frame_duration * fs);
	int Ns = round(0.001 * frame_shift * fs);

	vector <vector <float> > frames;
	vec2frames(speech, Nw, Ns, frames);
	return 0;
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

	vector <vector <float>> MFCCs_train;

	if (get_feature_vector(speech, header.sampleRate, MFCCs_train))
	{
		cout << "error get_feature_vector" << endl;
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
	fclose(file);

	_getch();
	return 0;
}