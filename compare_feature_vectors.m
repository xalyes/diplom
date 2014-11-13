function dist = compare_feature_vectors(MFCCs1, MFCCs2) %#codegen
% ������� compare_feature_vectors ���������� ������� ��������� MFCCs1 � MFCCs2
% ���������� dist - ��������� �������� ���������� ����� ��������� ���������
N1 = size(MFCCs1, 2);
N2 = size(MFCCs2, 2);
D = disteusq(MFCCs1', MFCCs2'); % ������� ���������� ����������
width = 0.5;
F = Inf(N1, N2); % �����
for i = 1:N1
for j = 1:N2
d = abs(i/N1 - j/N2);
if d <= width / 4
F(i, j) = 1;
else
if d <= width
F(i, j) = 1 + d * 5;
end;
end;
end;
end;
S = D.*F; % ������� ��������
dist = sqrt(sum(min(S, [], 2).^2) / size(S, 1));
%dist = sum(min(d, [], 2) / size(d, 1));
end