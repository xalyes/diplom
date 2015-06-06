#include "main.h";

using namespace std;
#define PI 3.1415926535
const double TwoPi = 6.283185307179586;

#define  NUMBER_IS_2_POW_K(x)   ((!((x)&((x)-1)))&&((x)>1))  // x is pow(2, k), k=1,2, ...
#define  FT_DIRECT        -1    // Direct transform.
#define  FT_INVERSE        1    // Inverse transform.

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

bool  FFT(float *Rdat, float *Idat, int N, int LogN, int Ft_Flag)
{
	// parameters error check:
	if ((Rdat == NULL) || (Idat == NULL))                  return false;
	if ((N > 16384) || (N < 1))                            return false;
	if (!NUMBER_IS_2_POW_K(N))                             return false;
	if ((LogN < 2) || (LogN > 14))                         return false;
	if ((Ft_Flag != FT_DIRECT) && (Ft_Flag != FT_INVERSE)) return false;

	register int  i, j, n, k, io, ie, in, nn;
	float         ru, iu, rtp, itp, rtq, itq, rw, iw, sr;

	static const float Rcoef[14] =
	{ -1.0000000000000000F, 0.0000000000000000F, 0.7071067811865475F,
	0.9238795325112867F, 0.9807852804032304F, 0.9951847266721969F,
	0.9987954562051724F, 0.9996988186962042F, 0.9999247018391445F,
	0.9999811752826011F, 0.9999952938095761F, 0.9999988234517018F,
	0.9999997058628822F, 0.9999999264657178F
	};
	static const float Icoef[14] =
	{ 0.0000000000000000F, -1.0000000000000000F, -0.7071067811865474F,
	-0.3826834323650897F, -0.1950903220161282F, -0.0980171403295606F,
	-0.0490676743274180F, -0.0245412285229122F, -0.0122715382857199F,
	-0.0061358846491544F, -0.0030679567629659F, -0.0015339801862847F,
	-0.0007669903187427F, -0.0003834951875714F
	};

	nn = N >> 1;
	ie = N;
	for (n = 1; n <= LogN; n++)
	{
		rw = Rcoef[LogN - n];
		iw = Icoef[LogN - n];
		if (Ft_Flag == FT_INVERSE) iw = -iw;
		in = ie >> 1;
		ru = 1.0F;
		iu = 0.0F;
		for (j = 0; j<in; j++)
		{
			for (i = j; i<N; i += ie)
			{
				io = i + in;
				rtp = Rdat[i] + Rdat[io];
				itp = Idat[i] + Idat[io];
				rtq = Rdat[i] - Rdat[io];
				itq = Idat[i] - Idat[io];
				Rdat[io] = rtq * ru - itq * iu;
				Idat[io] = itq * ru + rtq * iu;
				Rdat[i] = rtp;
				Idat[i] = itp;
			}

			sr = ru;
			ru = ru * rw - iu * iw;
			iu = iu * rw + sr * iw;
		}

		ie >>= 1;
	}

	for (j = i = 1; i<N; i++)
	{
		if (i < j)
		{
			io = i - 1;
			in = j - 1;
			rtp = Rdat[in];
			itp = Idat[in];
			Rdat[in] = Rdat[io];
			Idat[in] = Idat[io];
			Rdat[io] = rtp;
			Idat[io] = itp;
		}

		k = nn;

		while (k < j)
		{
			j = j - k;
			k >>= 1;
		}

		j = j + k;
	}

	if (Ft_Flag == FT_DIRECT) return true;

	rw = 1.0F / N;

	for (i = 0; i<N; i++)
	{
		Rdat[i] *= rw;
		Idat[i] *= rw;
	}

	return true;
}

void EraseColMatrinx(vector <vector <float>> &source, int N)
{
	for (int i = 0; i < source.size(); ++i)
	{
		source[i].erase(source[i].begin() + N);
	}
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
	if (A.front().size() != B.size())
		throw 1;
	int a = A.size();
	int b = B.size();
	int c = B.front().size();

	float **Am = new float*[a]; // ������ � �������
	for (int count = 0; count < a; count++)
		Am[count] = new float[b]; // �������
	float **Bm = new float*[b]; // ������ � �������
	for (int count = 0; count < b; count++)
		Bm[count] = new float[c]; // �������
	float **R = new float*[a]; // ������ � �������
	for (int count = 0; count < a; count++)
		R[count] = new float[c]; // �������
	vector <vector <float> > Rm(a, vector<float>(c, 0));

	for (int i = 0; i < a; ++i)
		for (int j = 0; j < b; ++j)
		{
			Am[i][j] = A[i][j];
		}
	for (int i = 0; i < b; ++i)
		for (int j = 0; j < c; ++j)
		{
			Bm[i][j] = B[i][j];
		}

	for (int i = 0; i < a; ++i)
		for (int j = 0; j < c; ++j)
		{
			R[i][j] = 0;
			for (int k = 0; k < b; ++k)
				R[i][j] += Am[i][k] * Bm[k][j];
		}
	for (int i = 0; i < a; ++i)
		for (int j = 0; j < c; ++j)
		{
			Rm[i][j] = R[i][j];
		}
	Res = Rm;

	for (int count = 0; count < a; count++)
		delete[] Am[count];
	for (int count = 0; count < b; count++)
		delete[] Bm[count];
	for (int count = 0; count < a; count++)
		delete[] R[count];
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

	vector <vector <float> > frames(inds.size(), vector<float>(indf.size()));
	for (int i = 0; i < inds.size(); ++i)
		for (int j = 0; j < indf.size(); ++j)
		{
			frames[i][j] = vec[inds[i] + indf[j] - 1];
		}
	
	vector<float> window_s(Nw);
	for (int i = 0; i < window_s.size(); ++i)
		window_s[i] = 0.54 - 0.46 * cos((2.0 * PI * i) / (window_s.size() - 1.0));
	vector<vector<float>> ham_matrix(Nw, vector<float>(Nw, 0));
	diag(window_s, ham_matrix);
	MulMatrix(ham_matrix, frames, resFrames);
}

void silence_removal(vector <vector <float>> &source)
{
	const int Nf = source.front().size();
	vector <float> E(Nf, 0);
	const int Nw = source.size();
	float E_mean = 0;
	for (int i = 0; i < Nf; ++i)
	{
		float sum = 0;
		for (int j = 0; j < Nw; ++j)
			sum += fabs(source[j][i])*fabs(source[j][i]);
		E[i] = sum * (1.0 / Nw);
		E_mean += E[i];
	}
	E_mean /= Nf;
	float T_E = E_mean / 10;
	vector<bool> voicedIndexes(Nf);
	for (int i = 0; i < E.size(); ++i)
		voicedIndexes[i] = ((E[i] > T_E) ? 1 : 0);
	for (int i = Nf - 1; i >= 0; --i)
		if (!voicedIndexes[i])
			EraseColMatrinx(source, i);

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

	//filter(25, speech);

	int Nw = round(0.001 * frame_duration * fs);
	int Ns = round(0.001 * frame_shift * fs);

	vector <vector <float> > frames;
	vec2frames(speech, Nw, Ns, frames);
	silence_removal(frames);
	const int Nfft = pow(2, (int)ceil(log((double)Nw) / log(2.0)));
	const int K = (Nfft / 2) + 1;
	int i;
	float *re = new float[Nfft];
	float *im = new float[Nfft];

	std::vector <std::vector <float>> MAG(Nfft, vector<float>(frames.front().size(), 0));
	for (int j = 0; j < frames.front().size(); ++j)
	{
		for (i = 0; i < Nfft; ++i)
		{
			re[i] = 0.0;
			im[i] = 0.0;
		}
		for (i = 0; i < 4410; ++i)
		{
			re[i] = frames[i][j];
		}

		FFT(re, im, Nfft, 13, FT_DIRECT);

		for (i = 0; i < Nfft; ++i)
		{
			MAG[i][j] = sqrt(re[i]*re[i] + im[i] * im[i]);
		}
	}
	delete[] re;
	delete[] im;
	return 0;
}

int _tmain(int argc, _TCHAR* argv[])
{
	FILE *file;
	errno_t err;
	err = fopen_s(&file, "train-01.wav", "rb");
	if (err)
	{
		printf_s("Failed open file, error %d", err);
		return 0;
	}

	WAVHEADER header;

	fread_s(&header, sizeof(WAVHEADER), sizeof(WAVHEADER), 1, file);
	vector<float> speech(header.subchunk2Size / 4);
	int j = 0; int i = 0;
	for (std::vector<float>::iterator it = speech.begin(); it != speech.end(); ++it)
	{
		float a, b;
		fread_s(&a, sizeof(a), sizeof(a), 1, file);
		//fread_s(&b, sizeof(b), sizeof(b), 1, file);
		//*it = bytesToFloat(a, b);
		*it = a;
		i++;
		if (i == 4410) j++;
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
	vector <vector <float>> MFCCs_train;

	if (get_feature_vector(speech, header.sampleRate, MFCCs_train))
	{
		cout << "error get_feature_vector" << endl;
	}
	
	fclose(file);

	_getch();
	return 0;
}