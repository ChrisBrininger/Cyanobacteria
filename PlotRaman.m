clc
clearvars

[baseName, folder] = uigetfile('*.txt','Select TxT File');
fullFileName = fullfile(folder, baseName);

fid = fopen(fullfile(fullFileName));

a = textscan(fid,'%f %f %*[^\n]', 'headerlines',50 );

x = a{1,1};

y = a{1,2};

plot(x,y)






