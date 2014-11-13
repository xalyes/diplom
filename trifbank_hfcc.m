function fbHfcc = trifbank_hfcc(Nhfcc, K, fmin, fmax, fs, ERBscaleFactor) %#codegen
% ������� fbHfcc = trifbank_hfcc(Nhfcc, K, fmin, fmax, fs, ERBscaleFactor)
%
% Mark Skowronski, March 26, 2003
Nfft = K * 2;
% ������������ ����������-����
a = 6.23e-6;
b = 93.39e-3;
c = 28.52;
a = a * ERBscaleFactor;
b = b * ERBscaleFactor;
c = c * ERBscaleFactor;
% ����������� ������� ������� �������
ahat = 0.5 * (1 / (700 + fmin));
bhat = 0.5 * (2 * 700 / (700 + fmin));
chat = 0.5 * (700 * (-1 + 700 / (700 + fmin)) - fmin);
bbar = (b - bhat) / (a - ahat);
cbar = (c - chat) / (a - ahat);
fcLow = 0.5 * (-bbar + sqrt(bbar^2 - 4 * cbar));
% ����������� ������� ������� �������
ahat = -0.5 * (1 / (700 + fmax));
bhat = -0.5 * (2 * 700 / (700 + fmax));
chat = 0.5 * (fmax - 700 * (-1 + 700 / (700 + fmax)));
bbar = (b - bhat) / (a - ahat);
cbar = (c - chat) / (a-ahat);
fcHigh = 0.5 * (-bbar + sqrt(bbar^2 - 4 * cbar));
% ��� �����
fcLowmel = 2595 * log10(1 + fcLow / 700);
fcHighmel = 2595 * log10(1 + fcHigh / 700);
melSpace = (fcHighmel - fcLowmel) / (Nhfcc - 1);
fcHfccmel = [0, fcLowmel + melSpace * (0:Nhfcc - 1), ...
2595 * log10(1 + (fs / 2) / 700)];
fcHfcc = 700 * (-1 + 10.^(fcHfccmel / 2595));
% ������������� ����� ��������
fbHfcc = zeros(length(fcHfcc) - 2, Nfft / 2);
fftFreqs = (0:Nfft / 2 - 1) / Nfft * fs;
% �������� ����� ��������
centerF = fcHfcc(2:length(fcHfcc) - 1);
centerFmel = 2595 * log10(1 + centerF / 700);
% ����������� ������ ������� �������
ERB = 6.23 * (centerF / 1000).^2 + 93.39 * (centerF / 1000) + 28.52;
% ��������������� ������ ��������
ERB = ERB * ERBscaleFactor;
% ������� ��������
lowerFmel = 2595 * log10(-2 / 1400 * ERB + 0.5 * sqrt((2 / 700 * ERB).^2 + ...
4 * 10.^(2 * centerFmel / 2595)));
upperFmel = 2 * centerFmel - lowerFmel;
% ������� ������ �������� � ����� ����
lowerF = 700 * (-1 + 10.^(lowerFmel / 2595));
upperF = 700 * (-1 + 10.^(upperFmel / 2595));
for chan = 1:length(centerF),
fbHfcc(chan, :) = (fftFreqs > lowerF(chan) & fftFreqs <= centerF(chan)).* ...
(fftFreqs - lowerF(chan)) / (centerF(chan) - lowerF(chan)) + ...
(fftFreqs > centerF(chan) & fftFreqs < upperF(chan)).*...
(upperF(chan) - fftFreqs) / (upperF(chan) - centerF(chan));
end
end