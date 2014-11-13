function dist = compare_feature_vectors(MFCCs1, MFCCs2) %#codegen
% Функция compare_feature_vectors сравнивает вектора признаков MFCCs1 и MFCCs2
% возвращает dist - численное значение расстояния между векторами признаков
N1 = size(MFCCs1, 2);
N2 = size(MFCCs2, 2);
D = disteusq(MFCCs1', MFCCs2'); % матрица евклидовых расстояний
width = 0.5;
F = Inf(N1, N2); % маска
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
S = D.*F; % матрица схожести
dist = sqrt(sum(min(S, [], 2).^2) / size(S, 1));
%dist = sum(min(d, [], 2) / size(d, 1));
end