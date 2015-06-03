#include "main.h";

using namespace std;
#define PI 3.1415926535
const double TwoPi = 6.283185307179586;

#define  NUMBER_IS_2_POW_K(x)   ((!((x)&((x)-1)))&&((x)>1))  // x is pow(2, k), k=1,2, ...
#define  FT_DIRECT        -1    // Direct transform.
#define  FT_INVERSE        1    // Inverse transform.

// Структура, описывающая заголовок WAV файла.
struct WAVHEADER
{
	// WAV-формат начинается с RIFF-заголовка:

	// Содержит символы "RIFF" в ASCII кодировке
	// (0x52494646 в big-endian представлении)
	char chunkId[4];

	// 36 + subchunk2Size, или более точно:
	// 4 + (8 + subchunk1Size) + (8 + subchunk2Size)
	// Это оставшийся размер цепочки, начиная с этой позиции.
	// Иначе говоря, это размер файла - 8, то есть,
	// исключены поля chunkId и chunkSize.
	unsigned long chunkSize;

	// Содержит символы "WAVE"
	// (0x57415645 в big-endian представлении)
	char format[4];

	// Формат "WAVE" состоит из двух подцепочек: "fmt " и "data":
	// Подцепочка "fmt " описывает формат звуковых данных:

	// Содержит символы "fmt "
	// (0x666d7420 в big-endian представлении)
	char subchunk1Id[4];

	// 16 для формата PCM.
	// Это оставшийся размер подцепочки, начиная с этой позиции.
	unsigned long subchunk1Size;

	// Аудио формат, полный список можно получить здесь http://audiocoding.ru/wav_formats.txt
	// Для PCM = 1 (то есть, Линейное квантование).
	// Значения, отличающиеся от 1, обозначают некоторый формат сжатия.
	unsigned short audioFormat;

	// Количество каналов. Моно = 1, Стерео = 2 и т.д.
	unsigned short numChannels;

	// Частота дискретизации. 8000 Гц, 44100 Гц и т.д.
	unsigned long sampleRate;

	// sampleRate * numChannels * bitsPerSample/8
	unsigned long byteRate;

	// numChannels * bitsPerSample/8
	// Количество байт для одного сэмпла, включая все каналы.
	unsigned short blockAlign;

	// Так называемая "глубиная" или точность звучания. 8 бит, 16 бит и т.д.
	unsigned short bitsPerSample;

	// Подцепочка "data" содержит аудио-данные и их размер.

	// Содержит символы "data"
	// (0x64617461 в big-endian представлении)
	char subchunk2Id[4];

	// numSamples * numChannels * bitsPerSample/8
	// Количество байт в области данных.
	unsigned long subchunk2Size;

	// Далее следуют непосредственно Wav данные.
};

float bytesToFloat(char firstByte, char secondByte)
{
	unsigned short int s = ((unsigned short int)secondByte << 8) + firstByte;
	float res;
	if (s >= 32768)
		res = (s - 65536) / 32768.0;
	else
		res = s / 32768.0;
	return floor(res * 100000) / 100000;
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
	vector <vector <float> > R(a, vector<float>(c, 0));

	for (int i = 0; i < a; ++i)
		for (int j = 0; j < c; ++j)
		{
			R[i][j] = 0;
			for (int k = 0; k < b; ++k)
				R[i][j] += A[i][k] * B[k][j];
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

void silence_removal(vector <vector <float>> source)
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
	for (int i = 0; i < Nf; ++i)
		if (!voicedIndexes[i])
			source.erase(source.begin() + i);

}

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

bool get_feature_vector(vector<float> speech, const int fs, vector <vector <float>> &MFCCs_train)
{
	bool use_ceplifter = false;
	float alpha = 0.97;				// коэффициент предварительной обработки
	int frame_duration = 100;		// длительность фрейма в мс
	int frame_shift = 50;			// сдвиг фреймов в мс
	int M = 40;						// количество фильтров
	int L = 22;						// параметр выравнивания
	float ERBScaleFactor = 0.75;	// коэффициент ширины фильтров
	int fmin = 100;					// минимальная частота, Гц
	int fmax = 6250;				// максимальная частота, Гц
	int Ncc = 29;					// количество мел - кепстральных коэффициентов

	//filter(25, speech);

	int Nw = round(0.001 * frame_duration * fs);
	int Ns = round(0.001 * frame_shift * fs);

	vector <vector <float> > frames;
	vec2frames(speech, Nw, Ns, frames);
	silence_removal(frames);
	const int Nfft = 8192;
	const int K = (Nfft / 2) + 1;
	float fram[8192];
	float res[8192] = {0.0};
	for (int j = 0; j < frames.size(); ++j)
	{
		for (int i = 0; i < 4410; ++i)
		{
			fram[i] = frames[i][j];
		}
		FFT(fram, res, Nfft, 13, FT_DIRECT);
	}
	return 0;
}

int _tmain(int argc, _TCHAR* argv[])
{
	/*vector<vector<float>> matrix1(2, vector<float>(3, 0));
	matrix1[0] = { 1, 2, 3 };
	matrix1[1] = { 2, 5, 10};
	vector<vector<float>> matrix2(3, vector<float>(2, 0));
	matrix2[0] = { 3, 4 };
	matrix2[1] = { 6, 7 };
	matrix2[2] = { 2, 1 };
	vector<vector<float>> matrixres(2, vector<float>(2, 0));
	MulMatrix(matrix1, matrix2, matrixres);*/

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
	// Выводим полученные данные
	printf_s("Sample rate: %d\n", header.sampleRate);
	printf_s("Channels: %d\n", header.numChannels);
	printf_s("Bits per sample: %d\n", header.bitsPerSample);
	printf_s("sizeof(header): %d\n", sizeof(short int));

	// Посчитаем длительность воспроизведения в секундах
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