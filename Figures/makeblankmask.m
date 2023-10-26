clearvars
clc

file = 'nd2/2021/channelbf,cy5,cfp,rfp_seq0000_0000.nd2';
ROI = [400 1525 400 400];

reader = BioformatsImage(file);


%%

for iT = 1:reader.sizeT

    I = getPlane(reader, 1, 'Cy5', iT);
    I = I(ROI(2):(ROI(2) + ROI(4)), ROI(1):(ROI(1) + ROI(3)));

    if iT == 1
        mask = false(size(I));


        imwrite(I, 'cell.tif', 'Compression', 'none')   
        imwrite(mask, 'mask_cell.tif', 'Compression', 'none');

    else
        
        imwrite(mask, 'mask_cell.tif', 'Compression', 'none', 'writemode', 'append');
        imwrite(I, 'cell.tif', 'Compression', 'none', 'writemode', 'append')
    end

end