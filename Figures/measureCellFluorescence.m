clearvars
clc

file = 'nd2/2021/channelbf,cy5,cfp,rfp_seq0000_0000.nd2';

maskFN = 'mask_cell.tif';

crystalMaskFN = 'mask_crystal.tif';

ROI = [400 1525 400 400];

L= LAPLinker;

L.LinkedBy = 'PixelIdxList';
L.LinkCostMetric = 'pxintersect';
L.LinkScoreRange = [1 12];
L.TrackDivision = true;
L.DivisionParameter = 'PixelIdxList';
L.DivisionScoreMetric = 'pxintersect';
L.DivisionScoreRange = [1 12];

reader = BioformatsImage(file);

dataCrystal = zeros(1,reader.sizeT);

for iT = 1:reader.sizeT

    I = getPlane(reader, 1, 'Cy5', iT);
    I = I(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)));

    mask = imread(maskFN, iT);
    mask = mask == 0;

    data = regionprops(mask, I, 'MeanIntensity', 'Area', 'PixelIdxList');

    L = assignToTrack(L, iT, data);

    crystalMask = imread(crystalMaskFN, iT);
    crystalMask = crystalMask > 0;
    tmpData = regionprops(crystalMask, 'Area');

    dataCrystal(iT) = tmpData.Area;

end

ct1 = getTrack(L, 1);
ct2 = getTrack(L, 2);
ct3 = getTrack(L, 3);

ts = reader.getTimestamps(1, 1);

tA = ts([ct1.Frames, ct2.Frames])/(60*60);
intA = [ct1.MeanIntensity; ct2.MeanIntensity];

tB = ts([ct1.Frames, ct3.Frames])/(60*60);
intB = [ct1.MeanIntensity; ct3.MeanIntensity];

tC = ts(1:numel(dataCrystal))/3600;
area = [dataCrystal] * prod(reader.pxSize);

figure(1)
plot(tA, intA, tB, intB)
xlim([0 60])

figure(2)
plot(tC, area)


