function HFCCs = get_feature_vector(speech, fs) %#codegen
% ������� get_feature_vector ��������� ������ (�������) ���������
% ��� ������������� �������� ������� speech � �������� ������������� fs
% ���������� HFCCs - ������ ���-������������ �������������
% ���� ������������� ����������
use_ceplifter = false;
alpha = 0.97; % ����������� ��������������� ���������
frame_duration = 100; % ������������ ������ � ��
frame_shift = 50; % ����� ������� � ��
M = 40; % ���������� ��������
L = 22; % �������� ������������
ERBScaleFactor = 0.75; % ����������� ������ ��������
fmin = 100; % ����������� �������, ��
fmax = 6250; % ������������ �������, ��
Ncc = 29; % ���������� ���-������������ �������������
% ����� ����� ������������� ����������
speech = filter(1 - alpha, 1, speech);
Nw = round(1E-3 * frame_duration * fs); % ����� ������ (� �������)
Ns = round(1E-3 * frame_shift * fs); % ����� ������ (� �������)
frames = vec2frames(speech, Nw, Ns); % ��������� ������� �� ������
voiced_frames = silence_removal(frames); % �������� ������
Nfft = 2 ^ nextpow2(Nw);
K = Nfft / 2 + 1;
MAG = abs(fft(voiced_frames, Nfft, 1)); % ���������� ������� �������
H = trifbank_hfcc(M, K, fmin, fmax, fs, ERBScaleFactor); % ��������� ����� ���-��������
FBE = H * MAG(1:K, :); % ���������� ���-��������
DCT = dctm(Ncc + 1, M); % ���������� ������� ����������� ��������������
CC = DCT * log(FBE); % ���������� ����������� ��������������
if use_ceplifter
lifter = ceplifter(Ncc + 1, L); % ���������� ������� ������������
figure;
plot(lifter);
HFCCs = diag(lifter) * CC; % ������������ �������� ���-������������ �������������
else
HFCCs = CC;
end;
HFCCs = HFCCs(2:end, :); % ������������ �������� ������������
end
function y = dctm(N, M)
% ������� ��������� ������� ����������� ��������������
y = sqrt(2.0 / M) * cos(repmat((0:N-1).', 1, M).* repmat(pi * ((1:M) - 0.5) / M, N, 1));
end
function y = ceplifter(N, L)
% ������� ������������ �������� ���-������������ �������������
y = 1 + 0.5 * L * sin(pi * (0:N-1) / L);
end