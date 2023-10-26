clearvars
clc

file = 'D:\Projects\Publication\2023 Brininger\20210418_Ana_-N_-Buffer_+CaCO3/channelbf,cy5,cfp,rfp_seq0000_0000.nd2';
maskFN = 'C:\Users\Jian Tay\OneDrive\Documents\Publications\2023 Brininger\code\export\fig2_c\crystalMask.tif';

%file = 'Z:\Chris Brininger\Movies\20220411_Ana_33047_H_Calcite_Bubble_-N/channelbf,cy5,cfp,rfp_seq0000_0000.nd2';

reader = BioformatsImage(file);


%ROI = [400 1525 400 400];
%ROI = [1647 1356 400 400]; 45:5:95

ROI = [1350 1420 400 400];

Icy5max = 11000;
Icy5min = 43;

cm = plasma(65536);
cm(1,:) = 0;

for iT = 1:reader.sizeT

    I = getPlane(reader, 1, 'Cy5', iT, 'ROI', ROI);
    Ibf = getPlane(reader, 1, 1, iT, 'ROI', ROI);

    % mask = imread(maskFN);
    % mask = mask == 0;
    
    Icy5norm = normalizeImg(I, Icy5max, Icy5min);

    Ibfmax = double(max(Ibf(:)));
    Ibfmin = double(min(Ibf(:)));

    Ibfnorm = normalizeImg(Ibf,  Ibfmax * 0.9, Ibfmin);

    % Icy5norm_red = zeros(size(I));
    % Icy5norm_red(:) = cm(Icy5norm(:) + 1, 1);
    % 
    % Icy5norm_green = zeros(size(I));
    % Icy5norm_green(:) = cm(Icy5norm(:) + 1, 2);
    % 
    % Icy5norm_blue = zeros(size(I));
    % Icy5norm_blue(:) = cm(Icy5norm(:) + 1, 3);

    % Icy5norm_red = cm(I(:) + 1, 1);
    % Icy5norm_red = imresize(Icy5norm_red, size(I));
    % 
    % Icy5norm_green = cm(I(:) + 1, 2);
    % Icy5norm_green = imresize(Icy5norm_green, size(I));
    % 
    % Icy5norm_blue = cm(I(:) + 1, 3);
    % Icy5norm_blue = imresize(Icy5norm_blue, size(I));

    % Icy5norm_rgb = cat(3, Icy5norm_red, Icy5norm_green, Icy5norm_blue);
    %imshow(Icy5norm_rgb)

    Ibfnorm = Ibfnorm ./65535;
    Icy5norm = Icy5norm ./65535;

    alpha = 0.7;

    % Imerge = cat(3, (alpha * Ibfnorm) + ((1 - alpha) * Icy5norm_rgb(:, :, 1)),...
    %     (alpha * Ibfnorm) + ((1 - alpha) * Icy5norm_rgb(:, :, 2)), ...
    %     (alpha * Ibfnorm) + ((1 - alpha) * Icy5norm_rgb(:, :, 3)));

    Imerge = cat(3, (alpha * Ibfnorm) + 3 * ((1 - alpha) * Icy5norm),...
        (alpha * Ibfnorm), ...
        (alpha * Ibfnorm));

    Imerge(Imerge > 1) = 1;
    Imerge(Imerge < 0) = 0;

    % imshow(Imerge(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)), :))
    % return

    Imerge = imresize(Imerge, 2);

    if iT == 1
        imwrite(Imerge, 'noGrowth.tif', 'Compression', 'none')
        imwrite(false(size(Imerge, 1), size(Imerge, 2)), 'crystalMask.tif', 'Compression', 'none')
    else
        imwrite(Imerge, 'noGrowth.tif', 'Compression', 'none', 'writeMode', 'append')
        imwrite(false(size(Imerge, 1), size(Imerge, 2)), 'crystalMask.tif', 'Compression', 'none', 'writeMode', 'append')
    end
      
    % imshow(Imerge)
    
    % imwrite(Imerge, ['Imerge_', int2str(iT), '.tif'], ...
    %      'Compression', 'none')


    % 
    % imwrite(Icy5norm,...
    %     ['Icy5_', int2str(iT), '.tif'], ...
    %     'Compression', 'none')
    % imwrite(Ibfnorm, ...
    %     ['Ibf_', int2str(iT), '.tif'], ...
    %     'Compression', 'none')
    % imwrite(Imerge, ...
    %     ['Imerge_', int2str(iT), '.tif'], ...
    %     'Compression', 'none')

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