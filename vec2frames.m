function frames = vec2frames(vec, Nw, Ns) %#codegen
% VEC2FRAMES Splits signal into overlapped frames using indexing.
%
% B=vec2frames(A,M,N) creates a matrix B whose columns consist of
% segments of length M, taken at every N samples along input vector A.
%
% [B,R]=vec2frames(A,M,N,D,W,P) creates a matrix B whose columns
% or rows, as specified by D, consist of segments of length M, taken
% at every N samples along the input vector A and windowed using the
% analysis window specified by W. The division of A into frames is
% achieved using indexes returned in R as follows: B=A(R);
%
% Summary
%
% A is an input vector
%
% M is a frame length (in samples)
%
% N is a frame shift (in samples)
% Author: Kamil Wojcicki, UTD, July 2011
% Modified by Mikhail Zadorozhnyy, 2014
vec = vec(:);
L = length(vec); % длина входного вектора
M = floor((L - Nw) / Ns + 1); % количество фреймов
% выравнивание
% определение необходимости выравнивания
E = (L - ((M - 1) * Ns + Nw));
% если выравнивание необходимо
if (E > 0)
P = Nw - E;
% дополнение нулями
vec1 = [ vec; zeros(P, 1) ];
M = M + 1;
else
vec1 = vec;
end
vec = vec1;
indf = Ns * (0:(M - 1)); % индексы фреймов
inds = (1:Nw).'; % индексы сэмплов
indexes = indf(ones(Nw, 1), :) + inds(:, ones(1, M)); % совмещенные индексы
frames = vec(indexes);
window_s = window_h(Nw);
frames = diag(window_s) * frames; % применение оконной функции
end
function y = window_h(N)
% функция окна Хэмминга
y = 0.54 - 0.46 * cos(2 * pi * (0:N-1).' / (N - 1));
end