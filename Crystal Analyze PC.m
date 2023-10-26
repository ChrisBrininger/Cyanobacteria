%% CyTracker to take a mask and a movie and get the data

CT = CyTracker;
importOptions (CT);
process(CT)

%%

File = BioformatsImage ( fullfile('E:\Movies\20210418_Ana_-N_-Buffer_+CaCO3\channelbf,cy5,cfp,rfp_seq0000_0000.nd2'));

ThickenAmount = [10,20,40,80,160,320];

Null = getPlane(File,1,1,1);
Null = (Null == 0);

CrystalCount = numel(tracks.Tracks);
% CrystalCount = 13;
FrameCount = numel(tracks.Tracks(1).Frames);

for iThick = 1:numel(ThickenAmount);
    
    for iCrystal = 1:CrystalCount
        %       iCrystal = 1
        
        if numel(tracks.Tracks(iCrystal).Data.PixelIdxList) == FrameCount
            for iFrame = 1:FrameCount
                %             iFrame = 2
                
                Cy5 = getPlane(File,1,2,iFrame);
                Mask = Null;
                %
                %
                Mask(tracks.Tracks(iCrystal).Data.PixelIdxList{1,iFrame}(:)) = 1;
                
                Area = bwmorph(Mask,'thicken',ThickenAmount(iThick));
                Area = Area - Mask;
                
                if numel(Area) == numel(Null)
                    %       RP = regionprops(Area,{‘Centroid’,‘PixelIdxList’});
                    RP = find(Area==1);
                    
                    %       Int = sum(Cy5(RP.PixelIdxList(:)));
                    Int = sum(Cy5(RP(:)))/numel(RP);
                    Intensity.Thicken{1,iThick}.Crystal{1,iCrystal}(iFrame) = Int;
                    
                else
                    
                    Intensity.Thicken{1,iThick}.Crystal{1,iCrystal}(iFrame) = [];
                    
                end
            end
        end
    end
end

%%
ThickenAmount = [10,20,40,80,160,320];
FrameCount = numel(tracks.Tracks(1).Frames);

X = 1:FrameCount;
close all

for iThick = 1:numel(ThickenAmount);
    %     iThick = 1
    figure
    hold on
    for iCrystal = 1:13
        
        %        iCrystal = 1
        if numel(tracks.Tracks(iCrystal).Data.PixelIdxList) == FrameCount
            
            Y = Intensity.Thicken{1,iThick}.Crystal{1,iCrystal}(:);
            
            if iCrystal == 1
                
                plot(X,Y,'Color','#fef162')
                
            elseif iCrystal == 2
                
                plot(X,Y,'Color','#e3b04f')
                
            elseif iCrystal == 3
                
                plot(X,Y,'Color','#f16474')
                
            elseif iCrystal == 4
                
                plot(X,Y,'Color','#c4da56')
                
            elseif iCrystal == 5
                
                plot(X,Y,'Color','#4c66b0')
                
            elseif iCrystal == 6
                
                plot(X,Y,'Color','#f49bc2')
                
            elseif iCrystal == 7
                
                plot(X,Y,'Color','#72c594')
                
            elseif iCrystal == 8
                
                plot(X,Y,'Color','#afe2f8')
                
            elseif iCrystal == 9
                
                plot(X,Y,'Color','#e86d3b')
                
            elseif iCrystal == 10
                
                plot(X,Y,'Color','#ecb121')
                
            elseif iCrystal == 11
                
                plot(X,Y,'Color','#e99449')
                
            elseif iCrystal == 12
                
                plot(X,Y,'Color','#635ca8')
                
            elseif iCrystal == 13
                
                plot(X,Y,'Color','#81cfd0')
                
            end
        end
    end
    
%     legend('Crystal 1','Crystal 2','Crystal 3','Crystal 4','Crystal 5','Crystal 6','Crystal 7','Crystal 8','Crystal 9','Crystal 10','Crystal 11','Crystal 12','Crystal 13');
    
end


%%

File = BioformatsImage ( fullfile('E:\Movies\20210418_Ana_-N_-Buffer_+CaCO3\channelbf,cy5,cfp,rfp_seq0000_0000.nd2'));

Null = getPlane(File,1,1,1);
Null = (Null == 0);

CrystalCount = numel(tracks.Tracks);
% CrystalCount = 13;
FrameCount = numel(tracks.Tracks(1).Frames);


for iCrystal = 1:CrystalCount
    %       iCrystal = 1
    
    if numel(tracks.Tracks(iCrystal).Data.PixelIdxList) == FrameCount
        iFrame =FrameCount;
        %             iFrame = 2
        
        Mask = Null;
        %
        %
        Mask(tracks.Tracks(iCrystal).Data.PixelIdxList{1,iFrame}(:)) = 1;
        
        imshow(Mask)
        pause;
        
    end
end




