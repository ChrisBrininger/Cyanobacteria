%% Thicken
Info = imfinfo('\Users\Chris\Documents\Boulder\MatLab\HeteroCysts\2020.3.5_ana33047_minusn_0003_Crop_series1_cellMask pre thicken.tif')

for iFrame = 1:numel(Info)
%     If the mask is being imported as black cells 
%     on white you need the imcomplement

    currMask = imread(Info(iFrame).Filename, iFrame);
%     currMask = imcomplement(currMask);
    currMask = bwmorph(currMask,'thicken',2);
 
    if iFrame == 1
        imwrite(currMask,'2020.3.5_ana33047_minusn_0003_Crop_series1_cellMask.tif','Compression','none')
    else
        imwrite(currMask,'2020.3.5_ana33047_minusn_0003_Crop_series1_cellMask.tif','WriteMode','append','Compression','none')
    end
    
    
end
%% CyTracker to take a mask and a movie and get the data

CT = CyTracker;

importOptions (CT);

process(CT)


%% Load Data
B = filamentDataAnalyzer('\Users\Chris\Documents\Boulder\Filament_Testing\Test\2020.3.5_ana33047_minusn_0003_Crop_series1.mat');


% Fix Time Stamps
%  B = fixTimeStamps (B, 0.166666666666666667)

%% Distance to Heterocysts 

% Note, the uM to nearest heterocyst is working I believe
B = NearestHeterocyst(B)

%%
ListofHet=[1,2,3]
E=52
ii=1

  for Hets = (ListofHet)
                                                                    
      Distance(end+1) = B.Tracks(Hets).Data.FilamentPosition{E}- B.Tracks(ii).Data.FilamentPosition{E}
 
  end


%% Plot

hold off

for iCell = 1:numel(B.Tracks)
    
    avgInt = nan(1,numel(B.Tracks(iCell).Frames));
    
    time = nan(1,numel(B.Tracks(iCell).Frames));
    
    for iFrame = 1:numel(B.Tracks(iCell).Frames)
        
    avgInt(iFrame) = ([B.Tracks(iCell).Data.TotalIntCy5{iFrame}]/numel(B.Tracks(iCell).Data.RegisteredPxInd{iFrame}))/([B.Tracks(iCell).Data.TotalIntCy5{1}]/numel(B.Tracks(iCell).Data.RegisteredPxInd{1}));
    
    time(iFrame)= B.Tracks(iCell).Frames(iFrame)*B.MeanDeltaT;

    end
    

    if avgInt(end)<=.80
%         && avgInt(end)>=3000
        
    
        disp(iCell)
    end
    
    if B.Tracks(iCell).isHeterocyst == 1
    plot( time,avgInt,'r')
    else
       plot( time,avgInt,'b')  
    % plot( B.Tracks(iCell).Frames,avgInt)
    end
    hold on
    
end

hold Off

%% Tree plot

treeplot(B,4)

%% Plot Lineage

plotLineage(B,24,'plotfit')

%% Plot the distance to nearest heterocyst in cells vs frame position

%plot(B.Tracks(1).Data.CellsToHeterocysts, B.Tracks(1).Frames)
cellnum=12;
plot(B.Tracks(cellnum).Frames, B.Tracks(cellnum).Data.MicronsToHeterocysts)
%% Plot the distance to nearest heterocyst in uM vs frame position

plotstyle = 2;
beforehetswap = 1;

maxFil = [];
maxDist = [];

%finds the greatest possible filament position
for iCell = 1:numel(B.Tracks)
    
    maxFil(end+1) = max ( cell2mat( B.Tracks(iCell).Data.FilamentPosition));
    
end

maxFil = max(maxFil);


%finds the greatest possible distance from a heterocyst
for iCell = 1:numel(B.Tracks)
    
    maxDist(end+1) = max ( B.Tracks(iCell).Data.MicronsToHeterocysts);
    
end

maxDist = max(maxDist);

% creates the blank matrix and fills with nan
DistMat = nan (numel(B.FileMetadata.Timestamps), maxFil);

for iFrame = 1:numel(B.FileMetadata.Timestamps)
    for iCell = 1:numel(B.Tracks)
        E = find(B.Tracks(iCell).Frames == iFrame);
        if isempty (E)==0
            DistMat (iFrame,B.Tracks(iCell).Data.FilamentPosition{E})=B.Tracks(iCell).Data.MicronsToHeterocysts(E);
        end
    end
    
 
  end

  for iFrame = 1:numel(B.FileMetadata.Timestamps)
      
       if   sum(DistMat(iFrame,:),"omitnan") == 0
      
     DistMat(iFrame,DistMat(iFrame,:)==0)=1200;
       end
  end

  
x = 1:maxFil;
y = 1:numel(B.FileMetadata.Timestamps);

figure(1)

[X,Y] = meshgrid (x,y);

% surf(X,Y,DistMat)
% colormap(jet)
if plotstyle == 1;
    
    surfl(X,(.166666*Y),DistMat)
    % colormap(green)    % change color map
    shading interp    % interpolate colors across lines and faces
    
elseif plotstyle == 2;
    
    surf(X,(.166666*Y),DistMat,'EdgeColor','none')
    mycolormap = customcolormap([0 .3 .9 1], {'#FFFFFF','#b9b8b9','#2bbec0','#b57855'});
    colormap(mycolormap);
%     shading interp    % interpolate colors across lines and faces
    colorbar;
    
elseif plotstyle == 3;
    
    surf(X,(.166666*Y),DistMat,'EdgeColor','none')
    mycolormap = customcolormap([0 .3 .9 1], {'#000000','#b4c0c8','#ea3c3c','#29abe2'});
    colormap(mycolormap);
    
%     shading interp    % interpolate colors across lines and faces
    colorbar;
    
elseif plotstyle == 4;
    
    surf(X,(.166666*Y),DistMat)
    mycolormap = customcolormap([0 .3 .9 1], {'#ffffff','#ffffff','#ff0000','#000000'});
    colormap(mycolormap);
    shading interp    % interpolate colors across lines and faces
    colorbar;
    
elseif plotstyle == 5;
    
    %plot distance between heterocysts over time as violin plot
    
    Dist = NaN (1,numel(DistMat(:,1)));
    
    for iFrame = 1:numel(DistMat(:,1))
        
        if   sum(DistMat(iFrame,:),"omitnan") == 0
            
            DistMat(iFrame,DistMat(iFrame,:)==0)=1200;
            
        else
            
            p = DistMat(iFrame,:);
            DistHet = findpeaks(p);
            DistHet = 2*mean(DistHet);
            
            Dist(1,iFrame) = DistHet;
            
        end
    end
    
    Dist(isnan(Dist))=max(Dist);
    
    plot(y,Dist)
    
end
% mesh(X,Y,DistMat)

%% Peak finder for above graph, enter frame number below and run along with
% the above script

FrameNum = (80);

%Plot the trace of what you are looking at, this is optional
figure(2)
plot (DistMat(FrameNum,:))

% Find the average and standard deviation of peaks
PKS = findpeaks(DistMat(FrameNum,:));
M = mean (PKS);
S = std(PKS);

disp ('mean maximum distance to het')
disp (M)
disp (char(177))
disp (S)

PKSx = PKS * 2;
M = mean (PKSx);
S = std(PKSx);

disp ('mean maximum distance between hets')
disp (M)
disp (char(177))
disp (S)

%%
% Find the ratio of heterocysts to normal cells in a given frame

%enter what frame number you want, or the range
FrameNum = (60:80);

maxFil = [];
maxDist = [];

%finds the greatest possible filament position
for iCell = 1:numel(B.Tracks)
    
    maxFil(end+1) = max ( cell2mat( B.Tracks(iCell).Data.FilamentPosition));
    
end

maxFil = max(maxFil);


%finds the greatest possible distance from a heterocyst
for iCell = 1:numel(B.Tracks)
    
    maxDist(end+1) = max ( B.Tracks(iCell).Data.MicronsToHeterocysts);
    
end

maxDist = max(maxDist);

% creates the blank matrix and fills with nan
DistMat = nan (numel(B.FileMetadata.Timestamps), maxFil);

for iFrame = 1:numel(B.FileMetadata.Timestamps)
    for iCell = 1:numel(B.Tracks)
        E = find(B.Tracks(iCell).Frames == iFrame);
        if isempty (E)==0
            DistMat (iFrame,B.Tracks(iCell).Data.FilamentPosition{E})=B.Tracks(iCell).Data.MicronsToHeterocysts(E);
        end
    end
end

x = 1:maxFil;
y = 1:numel(B.FileMetadata.Timestamps);

Ratios = [];

for iFrame = FrameNum
   
    Vall = findpeaks(max(DistMat(iFrame,:))-DistMat(iFrame,:));
    
    cellNumHet = [];
    cellNumNonHet = [];
    
    for iCell = 1:numel(B.Tracks)
        E = find(B.Tracks(iCell).Frames == iFrame);
        if isempty (E)==0
            if B.Tracks(iCell).isHeterocyst == 1
               
            else 
                cellNumNonHet(end+1) = B.Tracks (iCell).ID;
            end
        end
    end
    
    Ratios(end+1) = numel(Vall)/(numel (cellNumNonHet)+numel(Vall));
    
end

mean (Ratios)
std (Ratios)
plot (.1666*(60:80),Ratios)

%% Fit to new growth rates

NonHetGrwRate = [];

HetGrwRate = [];

cellNumHet = [];

cellNumNonHet = [];

cellFillPosHet = [];

cellFillPosNonHet = [];

AvgIntCells = [];

CellFrameCount = [];

for iCell = 1:numel(B.Tracks)
    
    AvgIntOverVid = [];
    
    if B.Tracks(iCell).isHeterocyst == 1
        
        HetGrwRate(end+1) = B.Tracks(iCell).GrowthRate;
        
        cellNumHet(end+1) = B.Tracks(iCell).ID;
        
        cellFillPosHet(end+1) = B.Tracks(iCell).Data.FilamentPosition{end};
        
    else 
        
        NonHetGrwRate(end+1) = B.Tracks(iCell).GrowthRate;
        
        cellNumNonHet(end+1) = B.Tracks (iCell).ID;
        
        cellFillPosNonHet(end+1) = B.Tracks(iCell).Data.FilamentPosition{end};
        
        CellFrameCount (end+1) = numel(B.Tracks(iCell).Frames);
        
        for iFrame = 1:numel(B.Tracks(iCell).Frames)
        
            AvgIntOverVid(end+1)=[B.Tracks(iCell).Data.TotalIntCy5{iFrame}]/numel(B.Tracks(iCell).Data.RegisteredPxInd{iFrame});
        end
        
        AvgIntCells(end+1) = sum (AvgIntOverVid)/numel(B.Tracks(iCell).Frames);
        
    end
       
    
    
end

% % --- Plot the cell number on x and the number of those cells on y --- % %

% hold Off
% 
% plot(numel(HetGrwRate),HetGrwRate,'b.')
% 
% hold On
% 
% plot(numel(NonHetGrwRate),NonHetGrwRate,'g.')
% 
% hold Off



% --- Plot the cell number on x and the growth rate on y --- % %
% 
% hold Off
% 
% plot(cellNumHet,HetGrwRate,'b.')
% 
% hold On
% 
% plot(cellNumNonHet,NonHetGrwRate,'g.')
% 
% hold Off
% 


% --- Plot the position on x and the growth rate on y --- % %

% hold Off
% 
% plot(cellFillPosHet,HetGrwRate,'b.')
% 
% hold On
% 
% plot(cellFillPosNonHet,NonHetGrwRate,'r.')
% 
% hold Off

% % --- Plot the cell number on x and the average intensity on y --- % %

% plot(cellNumNonHet,AvgIntCells,'g.')

% % --- Plot the cell frames count on x and the growth rate on y --- % %

% plot(CellFrameCount,NonHetGrwRate,'b.')

% % --- Plot the frame count on x and the average intensity on y --- % %

% plot(CellFrameCount,AvgIntCells,'b.')

% % --- Plot the average intensity on x and the growth rate on y --- % %

plot(AvgIntCells,NonHetGrwRate,'b.')

% % -- 3d Plot

%% 

hold on;

for iCrystal = (1:13)

x = .5 * (tracks.Tracks(1).Frames (1,:));

y = cell2mat(tracks.Tracks(iCrystal).Data.Area (1,:));

plot (x,y)

end

hold off

%%

currMaskDir = ('E:\Backup\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\2020.5.15 Anabaena 1 cell to 2 hets movie\Mask without false cell\ana33047_nh4_to_minusn_part1_0002_crop_crop(50-500 every 5th)_series1_cellMask.tif');
currMask = imread(currMaskDir,1);
imwrite(currMask, ('E:\Backup\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\2020.5.15 Anabaena 1 cell to 2 hets movie\Mask without false cell\ana33047_nh4_to_minusn_part1_0002_crop_crop(50-500 every 5th)_series1_cellMask_firstframe.tif'), 'compression', 'none');
currMask = imread(currMaskDir,72);
imwrite(currMask, ('E:\Backup\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\2020.5.15 Anabaena 1 cell to 2 hets movie\Mask without false cell\ana33047_nh4_to_minusn_part1_0002_crop_crop(50-500 every 5th)_series1_cellMask_finalframe.tif'), 'compression', 'none');

%%
clearvars
%Input the folder and generic file name (usually ending at seq000)
Folder = '20220823_gssu-roGFP-rbcL_gssu-roGFP_3LD_12hrcycles';
File = 'channelbf,gfp,cy5,rfp,bfp_seq000';
%Input what multipoint to combine
Multipoint = '0';
%Input as a list the number of segments to combine
Segments = [0 1 2 3 4 5];
%Adjust channel scaling if needed (default:1)
Channelscaling = [1 1 1 1 1];

for seg = 1:length(Segments)
%%Select file
    %Converts segment number into a character format for Bioformats
    Seg_num = char(string(Segments(seg)));
    %Creates full file name to access
    nd2file = strcat('F:\Widefield\', Folder, '\', File, Seg_num,'_000', Multipoint,'.nd2');
    %Create a new BioformatsImage object
    bfr = BioformatsImage(nd2file);
    %Find frame count and number of channels
    FrameCount = bfr.sizeT;
    NumChannels = length(bfr.channelNames);
    
    for iChannel = 1:NumChannels
    %%Cycle through each channel in the file
        %Pulls out channel name to save with
        Channel = char(bfr.channelNames(iChannel));
        %Adjusts channel intensity to desired scaling
        Scaling = Channelscaling(iChannel);
        %Creates full file name to save to
        tiffile = strcat('F:\Widefield\', Folder, '\', File, Multipoint, Channel, '_combined.tif');
        
        for iFrame = 1:FrameCount
        %%Cycle through each frame of the movie
            %Pulls out the channel from each frame
            currFrame = getPlane(bfr,1,iChannel,iFrame);
            %Creates the file to save to initially and then appends
            if iFrame == 1 && seg == 1
                imwrite((currFrame*Scaling), tiffile, 'compression', 'none');
            else
                imwrite((currFrame*Scaling), tiffile, 'writeMode', 'append', 'compression', 'none');
            end
        end
    end
end