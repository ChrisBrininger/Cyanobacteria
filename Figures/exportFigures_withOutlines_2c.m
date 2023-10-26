clearvars
clc

file = 'D:\Projects\Publication\2023 Brininger\20210418_Ana_-N_-Buffer_+CaCO3/channelbf,cy5,cfp,rfp_seq0000_0000.nd2';

outputDir = 'C:\Users\Jian Tay\OneDrive\Documents\Publications\2023 Brininger\code\export\fig2_c';

crystalMaskFN = 'C:\Users\Jian Tay\OneDrive\Documents\Publications\2023 Brininger\code\export\fig2_c\crystalMask.tif';

reader = BioformatsImage(file);

ROI = [1350 1420 400 400];

Icy5max = 11000;
Icy5min = 43;

for iT = 1:20:reader.sizeT

    I = getPlane(reader, 1, 'Cy5', iT, 'ROI', ROI);
    Ibf = getPlane(reader, 1, 1, iT, 'ROI', ROI);

    Icy5norm = normalizeImg(I, Icy5max, Icy5min);

    Ibfmax = double(max(Ibf(:)));
    Ibfmin = double(min(Ibf(:)));
    
    Ibfnorm = normalizeImg(Ibf,  Ibfmax * 0.7, Ibfmin);

    Ibfnorm = Ibfnorm ./65535;
    Icy5norm = Icy5norm ./65535;

    alpha = 0.7;

    Imerge = cat(3, (alpha * Ibfnorm) + 3 * ((1 - alpha) * Icy5norm),...
        (alpha * Ibfnorm), ...
        (alpha * Ibfnorm));

    Imerge(Imerge > 1) = 1;
    Imerge(Imerge < 0) = 0;

    Imerge = imresize(Imerge, 2);
    
    imshow(Imerge)
    hold on

    % cellMask = imread(maskFN, iT);
    % cellMask = cellMask == 0;
    % cellMask = imresize(cellMask, 2, 'nearest');

    crystalMask = imread(crystalMaskFN, iT);
    crystalMask = crystalMask > 0;
    % crystalMask = imresize(crystalMask, 2, 'nearest');
    
    % boundaryCell = bwboundaries(cellMask);
    boundaryCrystal = bwboundaries(crystalMask);

    % for iC = 1:numel(boundaryCell)
    %     plot(boundaryCell{iC}(:, 2), boundaryCell{iC}(:, 1), 'g')
    % end

    for iCr = 1:numel(boundaryCrystal)
        plot(boundaryCrystal{iCr}(:, 2), boundaryCrystal{iCr}(:, 1), 'b')
    end

    saveas(gcf, [outputDir, '\outlines_', int2str(iT), '.svg'])
    hold off


    % imwrite(imresize(Icy5norm(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)), :), 2),...
    %     ['Icy5_', int2str(iT), '.tif'], ...
    %     'Compression', 'none')
    % imwrite(imresize(Ibfnorm(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)), :), 2), ...
    %     ['Ibf_', int2str(iT), '.tif'], ...
    %     'Compression', 'none')
    % imwrite(imresize(Imerge(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)), :), 2), ...
    %     ['Imerge_', int2str(iT), '.tif'], ...
    %     'Compression', 'none')

end

function Inorm = normalizeImg(img, Imax, Imin)

Inorm = double(img);
Inorm = (Inorm - Imin)/(Imax - Imin);

Inorm(Inorm > 1) = 1;
Inorm(Inorm < 0) = 0;

Inorm = floor(Inorm * 65535);

end