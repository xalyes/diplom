#include "main.h";

using namespace std;
#define PI 3.1415926535

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

void filter(const float k, vector<float> &speech)
{
	vector<float> resSpeech(speech.size());
	int size = speech.size();
	for (int i = k; i < size - k; ++i)
	{
		float sum = 0;
		for (int j = -k; j < k + 1; ++j)
		{
			sum = speech[i + j];
		}
		resSpeech[i] = (1 / k)*sum;
	}
	for (int i = 0; i < k; ++i)
		resSpeech[speech.size() - i-1] = speech[i];
	for (int i = 0; i < k; ++i)
		resSpeech[i] = speech[i];
	speech = resSpeech;
}

void MulMatrix(const vector <vector <float> > &A, const vector <vector <float> > &B, vector <vector <float> > &Res)
{
	if (A.size() != B.front().size())
		return;
	int a = A.size();
	int b = A.size();
	int c = B.size();
	vector <vector <float> > R(c, vector<float>(a, 0));

	for (int i = 0; i < c; ++i)
		for (int j = 0; j < a; ++j)
		{
			R[i][j] = 0;
			for (int k = 0; k < b; ++k)
				R[i][j] += A[k][j] * B[i][k];
		}
	Res = R;
}

void diag(const vector <float> &V, vector <vector <float> > &R)
{
	vector <vector <float> > A(V.size(), vector<float>(V.size(), 0));
	for (int i = 0; i < V.size(); ++i)
		A[i][i] = V[i];
	R = A;
}

void vec2frames(vector<float> &vec, int Nw, int Ns, vector <vector <float> > &resFrames)
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

	vector <vector <float> > frames(indf.size(), vector<float>(inds.size()));
	for (int i = 0; i < indf.size(); ++i)
		for (int j = 0; j < inds.size(); ++j)
		{
			frames[i][j] = vec[indf[i] + inds[j]];
		}
	
	vector<float> window_s(Nw);
	for (int i = 0; i < window_s.size(); ++i)
		window_s[i] = 0.54 - 0.46 * cos((2.0 * PI * i) / (window_s.size() - 1.0));
	vector<vector<float>> ham_matrix(Nw, vector<float>(Nw, 0));
	diag(window_s, ham_matrix);
	MulMatrix(ham_matrix, frames, resFrames);
}

void silence_removal(vector <vector <float>> source)
{
	const int Nf = source.size();
	vector <float> E(Nf, 0);
	const int Nw = source.front().size();
	float E_mean = 0;
	for (int i = 0; i < Nf; ++i)
	{
		float sum = 0;
		for (int j = 0; j < Nw; ++j)
			sum += fabs(source[i][j])*fabs(source[i][j]);
		E[i] = sum * (1.0 / Nw);
		E_mean += E[i];
	}
	E_mean /= Nf;
	float T_E = E_mean / 10;
	vector<bool> voicedIndexes(Nf);
	for (int i = 0; i < E.size(); ++i)
		voicedIndexes[i] = ((E[i] > T_E) ? 1 : 0);
	for (int i = 0; i < Nf; ++i)
		if (!voicedIndexes[i])
			source.erase(source.begin() + i);

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

	filter(25, speech);

	int Nw = round(0.001 * frame_duration * fs);
	int Ns = round(0.001 * frame_shift * fs);

	vector <vector <float> > frames;
	vec2frames(speech, Nw, Ns, frames);
	silence_removal(frames);
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