%% CyTracker to take a mask and a movie and get the data

CT = CyTracker;

importOptions (CT);

process(CT)

%%

File =  BioformatsImage ( fullfile('/Users/chris/Desktop/FilesForAreaCrystal/channelbf,cy5,cfp,rfp_seq0000_0000.nd2'));

ThickenAmount = 100;

Null = getPlane(File,1,1,1);

Null = (Null == 0);

CrystalCount = numel(tracks.Tracks);

FrameCount = numel(tracks.Tracks(iCrystal).Frames);

for iCrystal = 1:CrystalCount
%     iCrystal = 1
    for iFrame = 1:FrameCount
%         iFrame = 2
        
        Cy5 = getPlane(File,1,2,iFrame);
        Mask = Null;
%         
%            
            Mask(tracks.Tracks(iCrystal).Data.PixelIdxList{1,iFrame}(:)) = 1;
            
            Area = bwmorph(Mask,'thicken',ThickenAmount);
            
            Area = Area - Mask;
            
%             RP = regionprops(Area,{'Centroid','PixelIdxList'});
RP = find(Area==1);

            
%             Int = sum(Cy5(RP.PixelIdxList(:)));
            
            Int = sum(Cy5(RP(:)));
            
            Intensity.Crystal{1,iCrystal}(iFrame) = Int;

    end



end
