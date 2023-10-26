clearvars
clc

reader = BioformatsImage('nd2/fig1/ana33047_nh4_to_minusn_part1_0002_crop_crop(50-500 every 5th).nd2');

ROI = [700 680 300 300];

vid = VideoWriter('movie1.avi');
vid.FrameRate = 5;
open(vid)

totalpxshift = [0 0];

for iT = 1:72

    Imask = imread('nd2/fig1/ana33047_nh4_to_minusn_part1_0002_crop_crop(50-500 every 5th)_series1_cellMask_edited.tif', iT);

    I = getPlane(reader, 1, 'Cy5', iT);
    Ibf = getPlane(reader, 1, 1, iT);

    Imasktmp = imdilate(Imask, strel('disk', 3));
    Imasktmp = Imasktmp > 1;
    data = regionprops(Imasktmp, 'Centroid');

    ROI = [data(1).Centroid(1) - 150, data(1).Centroid(2) - 150, ...
        300 300];
    ROI = round(ROI);

    if iT == 1

        % Iref = I;

        Ibfmax = 45000;
        Ibfmin = 6000;

        %Icy5max = double(max(I(:)));
        Icy5max = 13000;
        Icy5min = double(min(I(:)));

    else

        % pxshift = xcorrreg(Iref, I);
        % totalpxshift = pxshift;
        % 
        % I = shiftimg(I, totalpxshift);
        % Ibf = shiftimg(Ibf, totalpxshift);
        % Imask = shiftimg(Imask, totalpxshift);
        % 
        % Iref = I;

    end

    Icy5norm = normalizeImg(I, Icy5max, Icy5min);
    Ibfnorm = normalizeImg(Ibf,  Ibfmax, Ibfmin);

    alpha = 0.6;

    Imerge = cat(3, (alpha * Ibfnorm) + ((1 - alpha) * Icy5norm),...
        (alpha * Ibfnorm), (alpha * Ibfnorm));

    perimMask = bwperim(Imask);
    perimMask = imdilate(perimMask, [1 1; 1 1]);


    Iblend = showoverlay(Imerge, perimMask, 'Color', [1 1 0]);
    
    Iblend(Iblend > 1) = 1;
    Iblend(Iblend < 0) = 0;

    if any(iT == [1 21 41 61])

            imwrite(imresize(Icy5norm(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)), :), 2),...
                ['Icy5_', int2str(iT), '.tif'], ...
                'Compression', 'none')
            imwrite(imresize(Ibfnorm(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)), :), 2), ...
                ['Ibf_', int2str(iT), '.tif'], ...
                'Compression', 'none')
            imwrite(imresize(Imerge(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)), :), 2), ...
                ['Imerge_', int2str(iT), '.tif'], ...
                'Compression', 'none')

    end

    Icrop = Iblend(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)), :);

    Icrop = imresize(Icrop, 2, 'nearest');
    
    writeVideo(vid, Icrop)

end
close(vid)

function Inorm = normalizeImg(img, Imax, Imin)

    Inorm = double(img);
    Inorm = (Inorm - Imin)/(Imax - Imin);

    Inorm(Inorm > 1) = 1;
    Inorm(Inorm < 0) = 0;

end

