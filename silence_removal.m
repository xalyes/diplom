function frames = silence_removal(source) %#codegen
% Функция удаления тишины
% Отбрасывает фреймы, не содержащие полезного сигнала
% на основе анализ энергии фреймов
Nf = size(source, 2);
E = zeros(Nf, 1);
Nw = size(source, 1);
for i = 1:Nf
E(i) = (1 / Nw) * sum(abs(source(:, i)).^2);
end
E_mean = mean(E);
T_E = E_mean / 10;
voicedIndexes = E >= T_E;
frames = source(:, voicedIndexes);
end