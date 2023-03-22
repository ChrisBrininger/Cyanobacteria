%%  outputs number of pixels for 1, 10, 20, 30, etc
clearvars
clc

File =  BioformatsImage ( fullfile('E:\Movies\20210418_Ana_-N_-Buffer_+CaCO3\channelbf,cy5,cfp,rfp_seq0000_0000.nd2'));

FinalFrame = 101;
Results.AverageArea = [];

roiCFP = getPlane(File,1,3,FinalFrame);
% roiCFP = roiCFP < 225;
imshow(roiCFP,[])
colormap(parula)
pause;
roi = drawrectangle('StripeColor','y');
rect = [roi.Vertices(1,1),roi.Vertices(1,2),(roi.Vertices(3,1)-roi.Vertices(1,1)),(roi.Vertices(2,2)-roi.Vertices(1,2))];

iFrame = 1

currCFP = getPlane(File,1,3,(iFrame));

currCFProi = imcrop(currCFP,rect);
currCFProibw = currCFProi > 0;

imshow(currCFProi,[])

colormap(parula)

f = imfreehand;
% Then draw the region using imfreehand tool
bw = createMask(f);
outI = bw.*currCFProibw;
numPixelsInsideRegion = sum(outI(:));
%     roi = drawpolygon ('StripeColor','y');


Results.AverageArea(iFrame) = numPixelsInsideRegion;


for iFrame = 1:(FinalFrame)
    
    iFrame
    
    currCFP = getPlane(File,1,3,iFrame);
    
    currCFProi = imcrop(currCFP,rect);
    currCFProibw = currCFProi > 0;
    
    imshow(currCFProi,[])
    
    colormap(parula)
    
    f = imfreehand;
    % Then draw the region using imfreehand tool
    bw = createMask(f);
    outI = bw.*currCFProibw;
    numPixelsInsideRegion = sum(outI(:));
    %     roi = drawpolygon ('StripeColor','y');
    
    
    Results.AverageArea(iFrame) = numPixelsInsideRegion;
    
    
    
    
end

%% Auto-Seg %%

clearvars
clc



% %          Options           % %

FinalFrame = 101;


% %

File =  BioformatsImage ( fullfile('E:\Movies\20210418_Ana_-N_-Buffer_+CaCO3\channelbf,cy5,cfp,rfp_seq0000_0000.nd2'));

iFrame = 10;
% for iFrame = 1:(FinalFrame)
% 
CFP = getPlane(File,1,3,iFrame);


Cells = getPlane(File,1,2,iFrame);


% 

% histogram (CFP)
CFP = CFP > 225;
CFP = imopen(CFP, strel ('disk',5));
Cells = imbinarize (Cells);
Cells = imdilate(Cells, strel ('disk',5));
% CFP = imbinarize (CFP);

imshow(CFP,[])
% imshow(Cells,[])

Crystals = (CFP - Cells);
Crystals = (Crystals > 0);

Crystals = imopen(Crystals, strel ('disk',5));
Crystals = imclose(Crystals, strel ('disk',20));
Crystals = imdilate(Crystals, strel ('disk',5));
% Crystals = imopen(Crystals, strel ('disk',20));
rpCrystals = regionprops(Crystals,{'Centroid','MajorAxisLength','MinorAxisLength','Orientation','Area'});
imshow(Crystals)


























% histogram (currCFP)
% 
% CFP = currCFP;
% 
% 
% histogram (CFP)
% 
% 
% currCFProi = imcrop(CFP,rect);
% currCFProi = -(currCFProi-1);
% currCFProi = imbinarize (currCFProi);
% 
% currCFProi = imopen(currCFProi, strel ('disk',3));
% currCFProi= imclose(currCFProi,strel ('disk',20));
% 
% imshow(currCFProi)
% pause(1)
% 
% end
% CFP = imbinarize(currCFP);

% imshow(CFP)
% histogram(currCFP)

% 
% BF = imbinarize(currBF);
% BF = -(BF-1);
% BF = imclose(BF, strel ('disk',4));
% BF = imopen(BF, strel ('disk',6));
% % BF = imbinarize(BF);
% 
% 
% % Cy5 = imbinarize(currCy5);
% Cy5 = currCy5 > 2500;
% Cy5 = imopen(Cy5, strel ('disk',3));
% Cy5 = imdilate(Cy5,strel ('disk',1));
% 
% Hets = imbinarize(BF - Cy5);
% bHets = imopen(Hets, strel ('disk',5));
% Hets = imerode(bHets,strel ('disk',1));
% 
% % WSH = bwdist(~Cy5W);
% % WSH = -WSH;
% % WSH = imhmin(WSH,2); %remove shallow minima from the mask
% % WSH = watershed(WSH);
% 
% 
% 
% Cy5W = imerode (currCy5,strel ('disk',2));
% Cy5W = imbinarize(Cy5W);
% Cy5W = imerode(Cy5W,strel ('disk',0));
% 
% WS = bwdist(~Cy5W);
% WS = WS .* imgaussfilt(double(currCy5),3);%Need to apply smoothing to the cy5 channel
% WS = -WS;
% WS = imhmin(WS,2); %remove shallow minima from the mask
% WS = watershed(WS);
% 
% Cy5W = imdilate(Cy5W, strel('disk',5));
% 
% 
% Cy5W(WS == 0) = false;
% Mask = imbinarize(imbinarize(Cy5W - bHets) + Hets);
% nMask = Mask;
% 
% celldata = regionprops(Mask,'Area','Centroid','PixelList');
% 
% Areas = cat(celldata.Area);
% sCell = [];
% lCell = [];
% 
% for iCell = 1:numel(Areas)
% 
%     if Areas(iCell) <= 500;
%     
%         sCell(end+1) = iCell;
%         
%     elseif Areas(iCell) >= 1200;
%         
%         lCell(end+1) = iCell;
%         
%     end
% end
% 
% 
% 
% 
% % Lets delete every small cell, I want to combine every small cell with its
% % nearest centroid
% 
% % for iCell = 1:numel(sCell)
% %     
% %     n = Mask(celldata(sCell(iCell).PixelList));
% %     nMask(n)t = false;
% %     
% %     
% %     
% %     
% %     
% % end
% 
% 
% % x = [269 466];
% % y = [494 423];
% % z = improfile(currCy5,x,y);, grid on;
% 
% % figure(1);
% % imshow(BF,[])
% % % imshowpair (currBF, BF)
% % imshow(WS,[])
% figure(2);
% 
% imshowpair(bwperim(Mask),currBF)
% % % plot(z)
% % 
% % figure(3);
% % imshow(Hets,[])
% 
% % figure(4);
% 
% % I want to find out how to have the software ask me which cells are too
% % small if they should be combined with a cell next to it
% 
% figure(3)
% 
% histogram(cat(celldata.Area))
% 
% 
