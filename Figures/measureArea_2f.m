clearvars
clc

file = 'D:\Projects\Publication\2023 Brininger\20210418_Ana_-N_-Buffer_+CaCO3/channelbf,cy5,cfp,rfp_seq0000_0000.nd2';
maskFN = 'C:\Users\Jian Tay\OneDrive\Documents\Publications\2023 Brininger\code\export\fig2f\crystalMask.tif';

ROI = [1120 320 400 400];

reader = BioformatsImage(file);

dataCrystal = zeros(1,reader.sizeT);

for iT = 1:reader.sizeT

    crystalMask = imread(maskFN, iT);
    crystalMask = crystalMask > 0;
    tmpData = regionprops(crystalMask, 'Area');

    dataCrystal(iT) = tmpData.Area;

end
ts = reader.getTimestamps(1, 1);

tC = ts(1:numel(dataCrystal))/3600;
area = [dataCrystal] * prod(reader.pxSize./2);

figure(2)
plot(tC, area)
ylim([0 4000])



