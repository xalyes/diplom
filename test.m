% ������ ������������ �������� � �����������
%% ��������
disp('�������� �������');
train_path = 'train';
count = 2;
db_path = 'fv';
for i = 1:count
fname = sprintf('%s\\train-%02d.wav', train_path, i);
[speech, fs] = audioread(fname);
MFCCs_train = get_feature_vector(speech, fs);
fname1 = sprintf('%s\\%d.csv', db_path, i);
dlmwrite(fname1, MFCCs_train);
end;
db_path = 'fv';
train_count = 2;
threshold = 0;
% ����������� ����������� ��������
for i = 1:train_count-1
fname1 = sprintf('%s\\%d.csv', db_path, i);
MFCCs1 = dlmread(fname1);
for j = 2:train_count
if i == j
continue;
end;
fname2 = sprintf('%s\\%d.csv', db_path, j);
MFCCs2 = dlmread(fname2);
dist = compare_feature_vectors(MFCCs1, MFCCs2);
disp(['���������� ����� ��������� ��������� ', num2str(i), ' � ', num2str(j), ...
' ���������� ', num2str(dist)]);
threshold = threshold + dist;
end;
end;
c1 = 1.1;
threshold = threshold / count * 1.33;
threshold2 = threshold * c1;
disp(['��������� ���������� ������������� ', num2str(threshold)]);
disp(['��������� ���������� ���������� "�����������" ', num2str(threshold2)]);
%% �������� ����� ��������
disp('�������� ����� ��������');
verify_path = 'verify1';
wildcard = sprintf('%s/*.wav', verify_path);
files = dir(wildcard);
disp(['����� ��������: ', num2str(length(files))]);
mdist1 = 0;
for i = 1:length(files)
fname = sprintf('%s/%s', verify_path, files(i).name);
[speech, fs] = audioread(fname);
MFCCs = get_feature_vector(speech, fs);
distmin = inf;
for j = 1:count
fname1 = sprintf('%s\\%d.csv', db_path, j);
MFCCs_train = dlmread(fname1);
dist = compare_feature_vectors(MFCCs, MFCCs_train);
distmin = min(distmin, dist);
end;
if distmin < threshold
disp(['�������� ����� ', fname, ': ����������: ', num2str(distmin), ' ������: ������������']);
else
if distmin < threshold2
disp(['�������� ����� ', fname, ': ����������: ', num2str(distmin), ' ������: �����������']);
else
disp(['�������� ����� ', fname, ': ����������: ', num2str(distmin), ' ������: ����������']);
end;
end;
mdist1 = max(mdist1, distmin);
end;
%% �������� ����� �������� � �������������� ������
disp('�������� ����� ��������');
verify_path = 'verify2';
wildcard = sprintf('%s/*.wav', verify_path);
files = dir(wildcard);
disp(['����� ��������: ', num2str(length(files))]);
mdist2 = inf;
for i = 1:length(files)
fname = sprintf('%s/%s', verify_path, files(i).name);
[speech, fs] = audioread(fname);
MFCCs = get_feature_vector(speech, fs);
distmin = inf;
for j = 1:count
fname1 = sprintf('%s\\%d.csv', db_path, j);
MFCCs_train = dlmread(fname1);
dist = compare_feature_vectors(MFCCs, MFCCs_train);
distmin = min(distmin, dist);
end;
if distmin < threshold
disp(['�������� ����� ', fname, ': ����������: ', num2str(distmin), ' ������: ������������']);
else
if distmin < threshold2
disp(['�������� ����� ', fname, ': ����������: ', num2str(distmin), ' ������: �����������']);
else
disp(['�������� ����� ', fname, ': ����������: ', num2str(distmin), ' ������: ����������']);
end;
end;
mdist2 = min(mdist2, distmin);
end;
%% �������� �������������� ��������
disp('�������� �������������� ��������');
verify_path = 'verify3';
wildcard = sprintf('%s/*.wav', verify_path);
files = dir(wildcard);
disp(['����� ��������: ', num2str(length(files))]);
for i = 1:length(files)
fname = sprintf('%s/%s', verify_path, files(i).name);
[speech, fs] = audioread(fname);
MFCCs = get_feature_vector(speech, fs);
distmin = inf;
for j = 1:count
fname1 = sprintf('%s\\%d.csv', db_path, j);
MFCCs_train = dlmread(fname1);
dist = compare_feature_vectors(MFCCs, MFCCs_train);
distmin = min(distmin, dist);
end;
if distmin < threshold
disp(['�������� ����� ', fname, ': ����������: ', num2str(distmin), ' ������: ������������']);
else
if distmin < threshold2
disp(['�������� ����� ', fname, ': ����������: ', num2str(distmin), ' ������: �����������']);
else
disp(['�������� ����� ', fname, ': ����������: ', num2str(distmin), ' ������: ����������']);
end;
end;
end;
%% ������
dista = abs(mdist1 - mdist2);
distr = (mdist1 + mdist2) / 2;
kp = 1 + dista / (2 * distr);
disp(['������������ ���������� ��� ����� �������� ���������� ', num2str(mdist1)]);
disp(['����������� ���������� ��� ����� �������� ���������� ', num2str(mdist2)]);
disp(['������� ������������� ������ � ������������ ������ ���������� ���������� ', num2str(dista)]);
disp(['������������� ���������� ������������� ���������� ', num2str(distr)]);
disp(['������������� ����������� ������� "�����������" ���������� ', num2str(kp)]);