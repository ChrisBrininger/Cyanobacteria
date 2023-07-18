%% Geit Velocity
clearvars
clc
% % Options % %
initialFrame = 1;
finalFrame = 3;
File =  BioformatsImage ( fullfile('E:\Movies\20210418_Ana_-N_-Buffer_+CaCO3\channelbf,cy5,cfp,rfp_seq0000_0000.nd2'));
Time = 5; %(seconds)
pixelSize = 0.324; %(Microns)
% % Script % %
pFrame = (finalFrame - initialFrame) + 1;
Points = NaN(pFrame,2);
pFrame = pFrame - 1;
Velocity = NaN(pFrame,1);
for iFrame = (initialFrame):(finalFrame)
    Cells = getPlane(File,1,2,iFrame);
    imshow (Cells,[])
%     roi = drawrectangle('StripeColor','y');
%     rect = [roi.Vertices(1,1),roi.Vertices(1,2),(roi.Vertices(3,1)-roi.Vertices(1,1)),(roi.Vertices(2,2)-roi.Vertices(1,2))];
%
%     currROI = imcrop(Cells,rect);
%
%     imshow(currROI,[]);
%
    f = drawpoint;
    Points(iFrame,1) = f.Position(1,1);
    Points(iFrame,2) = f.Position(1,2);
    if iFrame == initialFrame;
    else
        dx = Points(iFrame,1) - Points(iFrame-1,1);
        dy = Points(iFrame,2) - Points(iFrame-1,2);
        H = sqrt(dx.^2 + dy.^2);
        H = H * pixelSize;
        Velocity(iFrame-1,1) = H/Time;
    end
end