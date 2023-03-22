% % User Inputs % %

% % % Note on limitations, currently the tif file can max out in size, so
% the script will limit the tif video to 350 frames per movie. This can
% output up to 3 movies currently, for a total of 1050 frames. The input
% for tif movie to mask will find and account for up to 3 movies.

% % These inputs are required for overlaying the mask % %

%Write the original mask directory here.
currDir = '\Users\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\2020.5.15 Anabaena 1 cell to 2 hets movie';

%Write name of the original mask here.
MaskFileName = 'ana33047_nh4_to_minusn_part1_0002_crop_crop(every10th)_shortened_series1_cellMask.tif';

%Write the name of the ND2 file to overlay here
BFFileName = 'ana33047_nh4_to_minusn_part1_0002_crop_crop(every10th).nd2';

%Write the desired name for the overlay mask (no numbers) here. 
%Must be the same name for output and input functions
maskName = 'maskOverlay';

%If the mask starts some number of frames in, enter the difference from 1
%here
CountDiff = 0;

%Write the directory for the output here.
outDir = '\Users\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\2020.5.15 Anabaena 1 cell to 2 hets movie\Test';

% % These additional inputs are required for going from overlay % %
% %                          back to mask                       % %

%Write the directory for the edited masks folder here.
maskDir = '\Users\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\2020.5.15 Anabaena 1 cell to 2 hets movie\Export\';

%Write the name you want given to the final recombined mask here.
outName = 'EditedMask';



%% Old Mask to Tif Movie Overlay %%

%First, load the metadata from the mask and start getting data from
%brightfield

MaskInfo = imfinfo(fullfile(currDir, MaskFileName));
BFFile =  BioformatsImage ( fullfile(currDir, BFFileName));

% This and the following if loops will determine how many frames the
% original mask has, and will use either 1, 2, or 3 tif movies as an
% output, capped at 350 frames per movie.

if numel(MaskInfo) <= 350
    
    % notes on what is happening only exist in this for loop.
    
    for iFrame = 1:numel(MaskInfo)
    
        % Call the old mask at the correct frame and make it into the
        % correct format (uint16). Additionally, take the complement
        % because Matlab makes no sense and will sometimes call the mask in
        % reverse.
        currMask = imread(fullfile(currDir, MaskFileName), iFrame);
        currMask = im2uint16(currMask);
        currMask = imcomplement (currMask);

        % Call the brightfield at the correct frame and format
        currBF = getPlane(BFFile,1,1,iFrame + CountDiff);

        % Overlay the mask as green and the brightfield as red and blue
        currOverlay = cat(3, currBF, currMask * .5, currBF);

         %Write to a new tif stack
        if iFrame == 1
            imwrite(currOverlay, fullfile(outDir,maskName,'1.tif'), 'compression', 'none');
        else
            imwrite(currOverlay, fullfile(outDir,maskName,'1.tif'), 'writeMode', 'append', 'compression', 'none');
        end

    end
    
else
    
    if numel(MaskInfo) <= 700
        for iFrame = 1:350
    
            currMask = imread(fullfile(currDir, MaskFileName), iFrame);
            currMask = im2uint16(currMask);
            currMask = imcomplement (currMask);

            currBF = getPlane(BFFile,1,1,iFrame + CountDiff);

            currOverlay = cat(3, currBF, currMask * .5, currBF);

             %Write to a new tif stack
            if iFrame == 1
                imwrite(currOverlay, fullfile(outDir,maskName,'1.tif'), 'compression', 'none');
            else
                imwrite(currOverlay, fullfile(outDir,maskName,'1.tif'), 'writeMode', 'append', 'compression', 'none');
            end

        end

        for iFrame = 351:numel(MaskInfo)

            currMask = imread(fullfile(currDir, MaskFileName), iFrame);
            currMask = im2uint16(currMask);
            currMask = imcomplement (currMask);

            currBF = getPlane(BFFile,1,1,iFrame + CountDiff);

            currOverlay = cat(3, currBF, currMask * .5, currBF);

             %Write to a new tif stack
            if iFrame == 1
                imwrite(currOverlay, fullfile(outDir,maskName,'2.tif'), 'compression', 'none');
            else
                imwrite(currOverlay, fullfile(outDir,maskName,'2.tif'), 'writeMode', 'append', 'compression', 'none');
            end

        end
        
    else
        
        if numel(MaskInfo) <= 1050
        
                for iFrame = 1:350
    
                    currMask = imread(fullfile(currDir, MaskFileName), iFrame);
                    currMask = im2uint16(currMask);
                    currMask = imcomplement (currMask);

                    currBF = getPlane(BFFile,1,1,iFrame + CountDiff);

                    currOverlay = cat(3, currBF, currMask * .5, currBF);

                     %Write to a new tif stack
                    if iFrame == 1
                        imwrite(currOverlay, fullfile(outDir,maskName,'1.tif'), 'compression', 'none');
                    else
                        imwrite(currOverlay, fullfile(outDir,maskName,'1.tif'), 'writeMode', 'append', 'compression', 'none');
                    end

                end

                for iFrame = 351:700
                    currMask = imread(fullfile(currDir, MaskFileName), iFrame);
                    currMask = im2uint16(currMask);
                    currMask = imcomplement (currMask);

                    currBF = getPlane(BFFile,1,1,iFrame + CountDiff);

                    currOverlay = cat(3, currBF, currMask * .5, currBF);

                     %Write to a new tif stack
                    if iFrame == 1
                        imwrite(currOverlay, fullfile(outDir,maskName,'2.tif'), 'compression', 'none');
                    else
                        imwrite(currOverlay, fullfile(outDir,maskName,'2.tif'), 'writeMode', 'append', 'compression', 'none');
                    end

                end
                
                for iFrame = 701:numel(MaskInfo)

                    currMask = imread(fullfile(currDir, MaskFileName), iFrame);
                    currMask = im2uint16(currMask);
                    currMask = imcomplement (currMask);

                    currBF = getPlane(BFFile,1,1,iFrame + CountDiff);

                    currOverlay = cat(3, currBF, currMask * .5, currBF);

                     %Write to a new tif stack
                    if iFrame == 1
                        imwrite(currOverlay, fullfile(outDir,maskName,'3.tif'), 'compression', 'none');
                    else
                        imwrite(currOverlay, fullfile(outDir,maskName,'3.tif'), 'writeMode', 'append', 'compression', 'none');
                    end

                end
        
        else
            
            error('Mask Greater than 1050 Frames.')
            
        end
        
    end
    
end
%% Tif movie Overlay to New Mask %% 

% Determine the number of movies %
numMovies = ( isfile (fullfile(maskDir,[maskName,sprintf('%01.0f',1),'.tif']))...
    +isfile (fullfile(maskDir,[maskName,sprintf('%01.0f',2),'.tif']))...
    +isfile (fullfile(maskDir,[maskName,sprintf('%01.0f',3),'.tif']))...
    +isfile (fullfile(maskDir,[maskName,sprintf('%01.0f',4),'.tif'])))


if numMovies <= 3
    
    for iMovie = 1:numMovies
    
        iMask = fullfile(maskDir,[maskName,sprintf('%01.0f',iMovie),'.tif']);
        EdittedMaskInfo = imfinfo(iMask);

        for iFrame = 1:numel (EdittedMaskInfo)

            % Calls the current file
            currMaskDir = fullfile(maskDir,[maskName,sprintf('%01.0f',iMovie),'.tif']);
            currMask = imread(currMaskDir, iFrame);

            % Creates the new binary mask based only on green.
            outputMask = currMask(:,:,2) > 0;

            %Write to a new tif stack
            if iMovie == 1
                if iFrame == 1
                    imwrite(outputMask, fullfile(outDir,[outName,'FromTif','.tif']), 'compression', 'none');
                else
                    imwrite(outputMask, fullfile(outDir,[outName,'FromTif','.tif']), 'writeMode', 'append', 'compression', 'none');
                end

            else

                imwrite(outputMask, fullfile(outDir,[outName,'FromTif','.tif']), 'writeMode', 'append', 'compression', 'none');

            end

        end


    end
    
else
    
   error('Currently can only handle a number of movies equal to 1, 2, or 3.')
   
end

%% Old Mask to PNG individual files Overlay %%

%First, load the metadata from the mask and start getting data from
%brightfield
MaskInfo = imfinfo(fullfile(currDir, MaskFileName));
BFFile =  BioformatsImage ( fullfile(currDir, BFFileName));
    
    for iFrame = 1:numel(MaskInfo)
    
    currMask = imread(fullfile(currDir, MaskFileName), iFrame);
    currMask = im2uint16(currMask);
    currMask = imcomplement (currMask);
    
    currBF = getPlane(BFFile,1,1,iFrame + CountDiff);
%     currBF = ((currBF - max(currBF(:))) / (max(currBF(:)) - min(currBF(:)))) * 65535;
    
    currOverlay = cat(3, currBF, currMask * .35, currBF);
    
     %Write to a new png named after the iFrame
     
     imwrite(currOverlay, [outDir,maskName,sprintf('%03.0f',iFrame),'.png'], 'compression', 'none');
   
    
end
    
%% PNG individual files to New Mask %%

%Takes the files, assuming the same number of files as there are frames 
% in the original mask, and takes the green pixels out, and turns them into
% a new binary mask named

% OldMaskInfo = imfinfo(fullfile(currDir, MaskFileName));

 currMaskDir = [maskDir,maskName,'001.png'];
    currMask = imread(currMaskDir);clc
    outputMask = currMask(:,:,2) > 0;
    imwrite(outputMask, [outName,'FromPNG','.tif'], 'compression', 'none');

    currMaskDir = [maskDir,maskName,'006.png'];
    currMask = imread(currMaskDir);clc
    outputMask = currMask(:,:,2) > 0;
    imwrite(outputMask, [outName,'FromPNG','.tif'], 'writeMode', 'append', 'compression', 'none');
    
for  iFrame = 1:35
    
    % Calls the current file
    currMaskDir = [maskDir,maskName,sprintf('%02.0f',iFrame),'1.png'];
    currMask = imread(currMaskDir);
    
    % Creates the new binary mask based only on green.
    outputMask = currMask(:,:,2) > 0;
    
    %Write to a new tif stack
%     if iFrame == 1
%         imwrite(outputMask, [outName,'FromPNG','.tif'], 'compression', 'none');
%     else
        imwrite(outputMask, [outName,'FromPNG','.tif'], 'writeMode', 'append', 'compression', 'none');
%     end
    
 % Calls the current file
    currMaskDir = [maskDir,maskName,sprintf('%02.0f',iFrame),'6.png'];
    currMask = imread(currMaskDir);
    
    % Creates the new binary mask based only on green.
    outputMask = currMask(:,:,2) > 0;
    
    %Write to a new tif stack
%     if iFrame == 1
%         imwrite(outputMask, [outName,'FromPNG','.tif'], 'compression', 'none');
%     else
        imwrite(outputMask, [outName,'FromPNG','.tif'], 'writeMode', 'append', 'compression', 'none');
%     end

end 