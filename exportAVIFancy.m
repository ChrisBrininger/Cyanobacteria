clearvars
clc
close all

%General Settings
% fnameND2 = '/Users/nicholashill/Desktop/seq0000_xy2_crop_middle 2 colonies.nd2';
% fnameMask = '/Users/nicholashill/Desktop/Matlab/Single Cell Analysis/ccmOp+ washouts/2nd exp N2P101/Middle 2 colonies/seq0000_xy2_crop_middle 2 colonies_series1_masks.tif';
fnameND2 = 'C:\Users\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\2020.3.5 Anabaena Original 80 frame movie and mask\Test\2020.3.5_ana33047_minusn_0003_Crop.nd2';
fnameMask = 'C:\Users\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\2020.3.5 Anabaena Original 80 frame movie and mask\2020.3.5_ana33047_minusn_0003_Crop_series1_cellMask.tif';
ROI = [300, 370, 200, 200]; %[XMIN YMIN WIDTH HEIGHT]
firstFrame = 1;
lastFrame = 80;

%Settings for BFandLabeledMasks version
channel = 'BF';

%Settings for FullImgandMaskandLabeledCBXandGraph
GFPLutRange = [58 820];
Cy5LutRange = [46 3200];
DoGSpotDiameter = 2.8;
SpotThreshold = 14;
MinSpotArea = 3;
ExtraFrames = 4; %So the movie does not loop immediately

%What version avi do you want?
version = 'BFandLabeledMasks';
%Version to choose from below
% % % % % % BFandLabeledMasks % % % % % % 
% % % % % % FullImgandMaskandLabeledCBXandGraph % % % % % %
% % % % % % FullImgWithArrows % % % % % %


%Video settings
outputFname = 'C:\Users\Chris\Documents\Boulder\Lab Work\Figures\movie';
vr_an = VideoWriter(outputFname, 'MPEG-4');
vr_an.Quality = 80;
vr_an.FrameRate = 4;

open(vr_an);

%Tracking settings
TL = LAPLinker;
TL.LinkedBy = 'PixelIdxList';
TL.LinkCostMetric = 'pxintersect';
TL.LinkScoreRange = [1, 12];
TL.MaxTrackAge = 2;
TL.TrackDivision = true;
TL.MinFramesBetweenDiv = 2;
TL.DivisionParameter = 'PixelIdxList';
TL.DivisionScoreMetric = 'pxintersect';
TL.DivisionScoreRange = [1, 12];
TL.Solver = 'lapjv';


%Here begins the code
bfr = BioformatsImage(fnameND2);

switch version
    
    case 'BFandLabeledMasks'
        
        for iT = firstFrame:lastFrame
            
            IBF = getPlane(bfr, 1, channel, iT, 'ROI', ROI);
            IBF = double(IBF)/65535;
            
            IBF = repmat(IBF, [1, 1, 3]);
            
            mask = imread(fnameMask, iT);
            mask = logical(mask);
            %     mask = mask(ROI(1):ROI(1)+ROI(4)-1, ROI(2):ROI(2)+ROI(3)-1);
            mask = mask(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
            
%             keyboard
            outlines = bwperim(mask);
            outlines = imdilate(outlines, ones(1));
            %Iout = showoverlay(double(I), outlines, 'normalize', true);
            Iout = uint8(mask) * uint8(255);
            
            cellData = regionprops(mask, {'Centroid', 'MajorAxisLength', 'PixelIdxList'});
            
            TL = assignToTrack(TL, iT, cellData);
            
            
            for iTrack = 1:TL.NumTracks
                
                currTrack = TL.tracks.Tracks(iTrack);
                
                if iT >= currTrack.Frames(1) && iT <= currTrack.Frames(end)
                    
                    trackCentroid = cat(2,currTrack.Data.Centroid);
                    
                    if isfield(currTrack, 'RegCentroid')
                        Iout = insertText(Iout, currTrack.Data.RegCentroid{end}, iTrack,...
                            'BoxOpacity', 0,'TextColor',[237, 85, 59]./255);
                    else
                        Iout = insertText(Iout, currTrack.Data.Centroid{end}, iTrack,...
                            'BoxOpacity', 0,'TextColor', 'blue', 'FontSize', 9, 'Font', 'LucidaBrightRegular', ...
                            'AnchorPoint', 'center');
                    end
                    
                    %             if iT > currTrack.FirstFrame
                    %                 Iout = insertShape(Iout, 'line', trackCentroid, 'color','white');
                    %             end
                    
                end
            end
            
            IBF = uint8(IBF * 255);
            
            writeVideo(vr_an, [IBF, Iout]);
            
            %     [A,map] = rgb2ind([IBF, Iout],256);
            %
            %
            %     if iT == 1
            %         imwrite(A,map,'movie.gif', 'Loopcount', inf, 'DelayTime', 0.14);
            %
            %     else
            %         imwrite(A,map, 'movie.gif', 'Writemode', 'append', 'DelayTime', 0.14);
            %     end
            
        end
        
        close(vr_an);
        

    case 'FullImgandMaskandLabeledCBXandGraph'
        
        for iT = firstFrame:lastFrame
            
            %Get GFP and Cy5 images
            GFP = getPlane(bfr, 1, 'GFP', iT, 'ROI', ROI); %* GFPIntensityMultiplier;
            Cy5 = getPlane(bfr, 1, 'Cy5', iT, 'ROI', ROI); %* Cy5IntensityMultiplier;
            GFP = double(GFP);
            Cy5 = double(Cy5);
            
            %Normalize images
            GFP = (GFP - GFPLutRange(1))/(GFPLutRange(2) - GFPLutRange(1));
            Cy5 = (Cy5 - Cy5LutRange(1))/(Cy5LutRange(2) - Cy5LutRange(1));
            
            %This will be leftmost image
            fullImg = cat(3, Cy5, GFP, GFP);
                        
            
            %Get mask
            mask = imread(fnameMask, iT);
            mask = logical(mask);
            mask = mask(ROI(2):ROI(2)+ROI(4)-1, ROI(1):ROI(1)+ROI(3)-1);
            
            %This will be middle image
            %maskFinalImg = uint8(mask) * uint8(255);
            maskFinalImg = repmat(mask, [1, 1, 3]);
            
            %create labeled image
            labeledImg = repmat(GFP, [1, 1, 3]);
            
            outlines = bwperim(mask);
            outlines = imdilate(outlines, ones(1));
            labeledImg = showoverlay(labeledImg, outlines, 'normalize', true, 'Color', [1 0 1]);
            
            
            
            %Add in the spot labels
            spotImg = getPlane(bfr, 1, 'GFP', iT, 'ROI', ROI);
            spotImg = double(spotImg);
            
            sigma1 = (1 / (1 + sqrt(2))) * DoGSpotDiameter;
            sigma2 = sqrt(2) * sigma1;
            
            g1 = imgaussfilt(spotImg, sigma1);
            g2 = imgaussfilt(spotImg, sigma2);
            
            dogImg = imcomplement(g2 - g1);
            
            %bgVal = mean(dogImg(:));
            
            [nCnts, xBins] = histcounts(dogImg(:));
            xBins = diff(xBins) + xBins(1:end-1);
            
            gf = fit(xBins', nCnts', 'gauss1');
            
            spotBg = gf.b1 + SpotThreshold .* gf.c1;
            
            %Segment the spots
            spotMask = dogImg > spotBg;
            spotMask(~mask) = false;
            
            spotMask = bwareaopen(spotMask, MinSpotArea);
            
            dd = -bwdist(~spotMask);
            
            LL = watershed(dd);
            
            %Final spot mask
            spotMask(LL == 0) = false;
            
            %From spot mask, outline each dot and overlay in image
            spotOutlines = bwperim(spotMask);
            spotOutlines = imdilate(spotOutlines, ones(1));
            labeledImg = showoverlay(labeledImg, spotOutlines, 'normalize', true, 'Color', [0 1 1]);
            
            
            %To assign track label
            cellData = regionprops(mask, {'Centroid', 'MajorAxisLength', 'PixelIdxList'});
            TL = assignToTrack(TL, iT, cellData);
            for iTrack = 1:TL.NumTracks
                
                currTrack = TL.tracks.Tracks(iTrack);
                
                if iT >= currTrack.Frames(1) && iT <= currTrack.Frames(end)
                    
                    trackCentroid = cat(2,currTrack.Data.Centroid);
                    
                    if isfield(currTrack, 'RegCentroid')
                        labeledImg = insertText(labeledImg, currTrack.Data.RegCentroid{end}, iTrack,...
                            'BoxOpacity', 0,'TextColor',[237, 85, 59]./255);
                    else
                        labelPos = currTrack.Data.Centroid{end};
                        labelPos = labelPos + 6;
                        labeledImg = insertText(labeledImg, labelPos, iTrack,...
                            'BoxOpacity', 0,'TextColor', 'magenta', 'FontSize', 18, 'Font', 'LucidaBrightRegular');
                    end
                    
                    %             if iT > currTrack.FirstFrame
                    %                 Iout = insertShape(Iout, 'line', trackCentroid, 'color','white');
                    %             end
                    
                end
            end
                    
            %Images require normalization to be written properly
            fullImg = fullImg / max(fullImg, [], 'all');
            fullImg(fullImg < 0) = 0;
            labeledImg = labeledImg / max(labeledImg, [], 'all');
            labeledImg(labeledImg < 0) = 0;
            
            
            %Make a spacer region between images
            spacer = ones(ROI(4), 10, 3);
            
            combinedImg = [fullImg, spacer, maskFinalImg, spacer, labeledImg];
            
            
            subplot(2, 1, 1);
            imshow(combinedImg)
           
            
            %Now, make the graph of cell length over time
            subplot(2, 1, 2);
            hold on
            xlim([1 21])
            xlabel('Frame')
            ylim([2 6])
            ylabel('Cell length (\mum)')
            
            colorVec = ['g', 'm', 'b', 'c', 'k', 'r', 'g'];
            allTexts = [];
            for iTrack = 1:numel(TL.tracks.Tracks)
                lengthVec = [TL.tracks.Tracks(iTrack).Data.MajorAxisLength{:}] * 0.065;
                timeVec = [TL.tracks.Tracks(iTrack).Frames];
                
                plot(timeVec, lengthVec, '.-', 'Color', colorVec(iTrack), 'LineWidth', 2, 'MarkerSize', 11)
                hText = text(timeVec(end) + 0.4, lengthVec(end), num2str(iTrack));
                num2add = numel(hText);
                allTexts(end + 1:end + num2add) = hText;
                
            end
            
            %Then, output the full image, including graph
            F = getframe(gcf);
            writeVideo(vr_an, F);
            if iT ~= lastFrame
                delete(allTexts)
            end
            

            
            
                
            %     [A,map] = rgb2ind([IBF, Iout],256);
            %
            %
            %     if iT == 1
            %         imwrite(A,map,'movie.gif', 'Loopcount', inf, 'DelayTime', 0.14);
            %
            %     else
            %         imwrite(A,map, 'movie.gif', 'Writemode', 'append', 'DelayTime', 0.14);
            %     end
            
        end
        
        for iFrame = 1:ExtraFrames
            writeVideo(vr_an, F);
        end
        
        close(vr_an);
        
        
        
        
    case 'FullImgWithArrows'
        
        
        
        
        
        



end %end of the switch statement



