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

	float **Am = new float*[a]; // строки в массиве
	for (int count = 0; count < a; count++)
		Am[count] = new float[b]; // столбцы
	float **Bm = new float*[b]; // строки в массиве
	for (int count = 0; count < b; count++)
		Bm[count] = new float[c]; // столбцы
	float **R = new float*[a]; // строки в массиве
	for (int count = 0; count < a; count++)
		R[count] = new float[c]; // столбцы
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

void trifbank_hfcc(unsigned Nhfcc, unsigned K, int fmin, int fmax, unsigned fs, float ERBscaleFactor, vector <vector <float>> &resFbHfcc)
{
	unsigned Nfft = K * 2;

	// Коэффициенты Глассберга - Мура
	float a = 6.23e-6;
	float b = 93.39e-3;
	float c = 28.52;
	a *= ERBscaleFactor;
	b *= ERBscaleFactor;
	c *= ERBscaleFactor;

	// Центральная частота первого фильтра
	float ahat = 0.5 * (1.0 / (700.0 + fmin));
	float bhat = 0.5 * (2.0 * 700.0 / (700.0 + fmin));
	float chat = 0.5 * (700.0 * (-1 + 700.0 / (700.0 + fmin)) - fmin);
	float bbar = (b - bhat) / (a - ahat);
	float cbar = (c - chat) / (a - ahat);
	float fcLow = 0.5 * (-bbar + sqrt(pow(bbar, 2) - 4.0 * cbar));

	// Центральная частота второго фильтра
	ahat = -0.5 * (1.0 / (700.0 + fmax));
	bhat = -0.5 * (2.0 * 700.0 / (700.0 + fmax));
	chat = 0.5 * (fmax - 700.0 * (-1 + 700.0 / (700.0 + fmax)));
	bbar = (b - bhat) / (a - ahat);
	cbar = (c - chat) / (a - ahat);
	float fcHigh = 0.5 * (-bbar + sqrt(pow(bbar, 2) - 4.0 * cbar));

	// Мел шкала
	float fcLowmel = 2595.0 * log10(1 + fcLow / 700.0);
	float fcHighmel = 2595.0 * log10(1 + fcHigh / 700.0);
	float melSpace = (fcHighmel - fcLowmel) / (Nhfcc - 1);
	vector <float> fcHfccmel(Nhfcc + 2);
	fcHfccmel[0] = 0;
	for (int i = 0; i < Nhfcc; ++i)
		fcHfccmel[i + 1] = fcLowmel + melSpace * i;
	fcHfccmel[Nhfcc+1] = 2595.0 * log10(1.0 + (fs / 2.0) / 700.0);
	vector <float> fcHfcc(fcHfccmel.size());
	for (int i = 0; i < fcHfcc.size(); ++i)
		fcHfcc[i] = 700.0 * (-1.0 + pow(10, fcHfccmel[i] / 2595.0));

	// Инициализация банка фильтров
	vector<vector<float>> fbHfcc(fcHfcc.size()-2, vector<float>(Nfft / 2, 0));
	vector<float> fftFreqs(Nfft / 2);
	for (int i = 0; i < Nfft / 2; ++i)
		fftFreqs[i] = ((i / (2.0 - 1.0)) / Nfft) * fs;

	// Создание банка фильтров
	vector <float> centerF(fcHfcc.size()-2);
	vector <float> centerFmel(fcHfcc.size() - 2);
	for (int i = 2; i < fcHfcc.size(); ++i)
	{
		centerF[i-2] = fcHfcc[i - 1];
		centerFmel[i-2] = 2595.0 * log10(1.0 + centerF[i-2] / 700.0);
	}

	// Определение ширины каждого фильтра и масштабирование ширины
	vector <float> ERB(centerF.size());
	for (int i = 0; i < ERB.size(); ++i)
	{
		ERB[i] = (6.23 * pow((centerF[i] / 1000.0), 2) + 93.39 * (centerF[i] / 1000.0) + 28.52)*ERBscaleFactor;
	}

	// Границы фильтров
	vector<float> lowerFmel(centerF.size());
	vector<float> upperFmel(centerF.size());
	for (int i = 0; i < lowerFmel.size(); ++i)
	{
		lowerFmel[i] = 2595.0 * log10((-2.0 / 1400.0) * ERB[i] + 0.5 * sqrt(pow(((2.0 / 700.0) * ERB[i]), 2) + 4.0 * pow(10, (2.0 * centerFmel[i] / 2595.0))));
		upperFmel[i] = 2.0 * centerFmel[i] - lowerFmel[i];
	}
	
	// Перевод границ фильтров в шкалу герц
	vector<float> lowerF(centerF.size());
	vector<float> upperF(centerF.size());
	for (int i = 0; i < lowerFmel.size(); ++i)
	{
		lowerF[i] = 700.0 * (-1.0 + pow(10.0, (lowerFmel[i] / 2595.0)));
		upperF[i] = 700.0 * (-1.0 + pow(10.0, (upperFmel[i] / 2595.0)));
	}

	vector<bool> freq(fbHfcc.front().size());
	vector<bool> freq2(fbHfcc.front().size());
	for (int chan = 0; chan < centerF.size(); ++chan)
	{
		for (int j = 0; j < fbHfcc.front().size(); ++j)
		{
			freq[j] = (fftFreqs[j] > lowerF[chan]) && (fftFreqs[j] <= centerF[chan]);
			freq2[j] = (fftFreqs[j] > centerF[chan]) && (fftFreqs[j] < upperF[chan]);
			fbHfcc[chan][j] = freq[j] * (fftFreqs[j] - lowerF[chan]) / (centerF[chan] - lowerF[chan]) \
				+ freq2[j] * (upperF[chan] - fftFreqs[j]) / (upperF[chan] - centerF[chan]);
		}
	}
	resFbHfcc = fbHfcc;
}

void dctm(int N, int M, vector <vector <float>> &DCT)
{
	vector <vector <float>> matrix1(N, vector<float>(M, 0));
	vector <vector <float>> matrix2(N, vector<float>(M, 0));
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
		{
			matrix1[i][j] = i;
			matrix2[i][j] = PI*(j+1.0-0.5) / M;
			DCT[i][j] = cos(matrix1[i][j] * matrix2[i][j]) * sqrt(2.0 / M);
		}
}

void ceplifter(int N, int L, vector <float> &y)
{
	for (int i = 0; i < N; ++i)
		y[i] = 1 + 0.5 * L * sin(PI * i / L);
}

void get_feature_vector(vector<float> speech, const int fs, vector <vector <float>> &MFCCs_train)
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
	const int Nfft = pow(2, (int)ceil(log2(Nw)));
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
	
	vector<vector<float>> H(M, vector<float>(K, 0));
	trifbank_hfcc(M, K, fmin, fmax, fs, ERBScaleFactor, H);

	vector<vector<float>> FBE(H.size(), vector<float>(MAG.front().size(), 0));

	MAG.erase(MAG.begin() + K, MAG.end());
	MulMatrix(H, MAG, FBE);

	vector<vector<float>> DCT(Ncc + 1, vector<float> (M, 0));
	dctm(Ncc + 1, M, DCT);

	vector<vector<float>> CC(DCT.size(), vector<float>(FBE.front().size(), 0));
	for (int i = 0; i < FBE.size(); ++i)
		for (int j = 0; j < FBE.front().size(); ++j)
			FBE[i][j] = log(FBE[i][j]);
	MulMatrix(DCT, FBE, CC);
	vector<vector<float>> HFCCs = CC;
	HFCCs.erase(HFCCs.begin());
	MFCCs_train = HFCCs;
}

void disteusq(vector <vector <float>> X, vector <vector <float>> Y, vector <vector <float>> &D)
{
	vector <vector<float>> x(X.front().size(), vector<float>(X.size()));
	vector <vector<float>> y(Y.front().size(), vector<float>(Y.size()));
	for (int i = 0; i < x.size(); ++i)
		for (int j = 0; j < x.front().size(); ++j)
			x[i][j] = X[j][i];
	for (int i = 0; i < y.size(); ++i)
		for (int j = 0; j < y.front().size(); ++j)
			y[i][j] = Y[j][i];
	int nx = x.size();
	int ny = y.size();
	int p = x.front().size();
	vector <vector<float>> d(nx, vector<float>(ny, 0));
	if (nx != ny)
	{
		vector <vector <vector<float>>> z(nx, vector<vector<float>>(ny, vector<float>(p)));
		vector <vector <vector<float>>> zy(nx, vector<vector<float>>(ny, vector<float>(p)));
		vector <vector<float>> d(nx, vector<float>(ny, 0));
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < p; ++j)
				for (int k = 0; k < ny; ++k)
				{
					z[i][k][j] = x[i][j];
				}
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < ny; ++j)
				for (int k = 0; k < p; ++k)
				{
					zy[i][j][k] = x[i][j];
					z[i][j][k] -= zy[i][j][k];
					d[i][j] += z[i][j][k] * z[i][j][k];
				}
		D = d;
	}
	else
	{
		for (int i = 0; i < nx; ++i)
			for (int j = 0; j < p; ++j)
			{
				d[i][j] = (x[i][j] - y[i][j]) * (x[i][j] - y[i][j]);
			}
			D = d;
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	FILE *file;
	errno_t err;
	long sampleRate1, sampleRate2;
	err = fopen_s(&file, "train-01.wav", "rb");
	if (err)
	{
		printf_s("Failed open file, error %d", err);
		return 0;
	}

	WAVHEADER header;

	fread_s(&header, sizeof(WAVHEADER), sizeof(WAVHEADER), 1, file);
	vector<float> speech1(header.subchunk2Size / 4);
	int j = 0; int i = 0;
	for (std::vector<float>::iterator it = speech1.begin(); it != speech1.end(); ++it)
	{
		float a, b;
		fread_s(&a, sizeof(a), sizeof(a), 1, file);
		*it = a;
		i++;
		if (i == 4410) j++;
	}
	sampleRate1 = header.sampleRate;
	err = fopen_s(&file, "train-02.wav", "rb");
	if (err)
	{
		printf_s("Failed open file, error %d", err);
		return 0;
	}

	fread_s(&header, sizeof(WAVHEADER), sizeof(WAVHEADER), 1, file);
	vector<float> speech2(header.subchunk2Size / 4);
	j = 0; i = 0;
	for (std::vector<float>::iterator it = speech2.begin(); it != speech2.end(); ++it)
	{
		float a, b;
		fread_s(&a, sizeof(a), sizeof(a), 1, file);
		*it = a;
		i++;
		if (i == 4410) j++;
	}
	sampleRate2 = header.sampleRate;
	vector <vector <float>> MFCCs_train1;
	vector <vector <float>> MFCCs_train2; // (28, vector<float>(2));

	get_feature_vector(speech1, sampleRate1, MFCCs_train1);
	get_feature_vector(speech2, sampleRate2, MFCCs_train2);

	
	vector <vector <float>> D;
	int N1 = MFCCs_train1.front().size();
	int N2 = MFCCs_train2.front().size();
	disteusq(MFCCs_train1, MFCCs_train2, D);
	float width = 0.5;
	float d;
	vector <vector <float>> F(N1, vector<float>(N2, FLT_MAX));
	for (i = 0; i < N1; ++i)
		for (j = 0; j < N2; ++j)
		{
			d = fabs((float) i/N1 - j/N2);
			if (d <= width / 4)
				F[i][j] = 1.0;
			else
				if (d <= width)
					F[i][j] = 1 + d * 5;
		}
	vector <vector <float>> S(N1, vector<float>(N2, 0));
	for (i = 0; i < N1; ++i)
		for (j = 0; j < N2; ++j)
			S[i][j] = D[i][j] * F[i][j];
	vector<float> minS(N1, FLT_MAX);
	float sum = 0;
	for (i = 0; i < N1; ++i)
	{
		for (j = 0; j < N2; ++j)
			if (S[i][j] < minS[i])
				minS[i] = S[i][j];
		sum += minS[i] * minS[i];
	}
	float dist = sqrt(sum / S.size());

	fclose(file);

	_getch();
	return 0;
}