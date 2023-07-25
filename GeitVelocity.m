%% Geit Velocity
%% Auto-Seg %%

clearvars
clc

% % Optons % %

initFrame = 1;
fFrame = 5;     %# of frames
File = BioformatsImage(fullfile('E:\Movies\20210418_Ana_-N_-Buffer_+CaCO3\channelbf,cy5,cfp,rfp_seq0000_0000.nd2')); %File location here
timeBtwnPhoto = 5;  %in seconds (get this from metadata!!!)
pixelSize = 0.324; %microns (get this from metadata!!!)
%remember, this is now in microns/sec!

% % Script % %

pFrame = fFrame - initFrame + 1;
points = NaN(pFrame, 2); %creates a pFrame x 2 vector (really an array but whos asking)
velocity = NaN(1,pFrame-1);
velIter = 1;
for iFrame = (initFrame):(fFrame)
    
    close all
    
    disp(iFrame)
    Cells = getPlane(File,1,1,iFrame);
    imshow(Cells,[])
    colormap('parula');
    f = drawpoint;
    %need to consider pre-allocation, or switching to a lang that is more
    %flexible...
    points(iFrame,1) = f.Position(1,1);
    points(iFrame,2) = f.Position(1,2);
    if iFrame ~= initFrame
        %Calculate velocity
        %find dist btwn two points using pythag
        dx = points(iFrame,1) - points(iFrame-1,1);
        dy = points(iFrame,2) - points(iFrame-1,2);
        dist = sqrt((dx.^2) + (dy.^2));
        dist = dist * pixelSize;
        velocity(velIter) = dist/timeBtwnPhoto;
        velIter = velIter + 1;
    end
end

