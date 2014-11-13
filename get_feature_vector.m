function HFCCs = get_feature_vector(speech, fs) %#codegen
% Функция get_feature_vector вычисляет вектор (матрицу) признаков
% для оцифрованного речевого сигнала speech с частотой дискретизации fs
% Возвращает HFCCs - марицу мел-кепстральных коэффициентов
% Блок настраиваемых параметров
use_ceplifter = false;
alpha = 0.97; % коэффициент предварительной обработки
frame_duration = 100; % длительность фрейма в мс
frame_shift = 50; % сдвиг фреймов в мс
M = 40; % количество фильтров
L = 22; % параметр выравнивания
ERBScaleFactor = 0.75; % коэффициент ширины фильтров
fmin = 100; % минимальная частота, Гц
fmax = 6250; % максимальная частота, Гц
Ncc = 29; % количество мел-кепстральных коэффициентов
% Конец блока настраиваемых параметров
speech = filter(1 - alpha, 1, speech);
Nw = round(1E-3 * frame_duration * fs); % длина фрейма (в сэмплах)
Ns = round(1E-3 * frame_shift * fs); % сдвиг фрейма (в сэмплах)
frames = vec2frames(speech, Nw, Ns); % Разбиение сигнала на фреймы
voiced_frames = silence_removal(frames); % Удаление тишины
Nfft = 2 ^ nextpow2(Nw);
K = Nfft / 2 + 1;
MAG = abs(fft(voiced_frames, Nfft, 1)); % Вычисление спектра сигнала
H = trifbank_hfcc(M, K, fmin, fmax, fs, ERBScaleFactor); % Получение банка мел-фильтров
FBE = H * MAG(1:K, :); % Применение мел-фильтров
DCT = dctm(Ncc + 1, M); % Вычисление матрицы косинусного преобразования
CC = DCT * log(FBE); % Применение косинусного преобразования
if use_ceplifter
lifter = ceplifter(Ncc + 1, L); % Вычисление матрицы выравнивания
figure;
plot(lifter);
HFCCs = diag(lifter) * CC; % Выравнивание значений мел-кепстральных коэффициентов
else
HFCCs = CC;
end;
HFCCs = HFCCs(2:end, :); % Отбрасывание старшего коэффициента
end
function y = dctm(N, M)
% Функция постоения матрицы косинусного преобразования
y = sqrt(2.0 / M) * cos(repmat((0:N-1).', 1, M).* repmat(pi * ((1:M) - 0.5) / M, N, 1));
end
function y = ceplifter(N, L)
% Функция выравнивания значений мел-кепстральных коэффициентов
y = 1 + 0.5 * L * sin(pi * (0:N-1) / L);
end