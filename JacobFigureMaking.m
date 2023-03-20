% % The goal of this script is to overlay the GFP channel onto the
% brightfield channel and then the outline of the mask, followed by the
% cell number from analysis.



    % % Settings % %

%Write the original mask directory here.
currDir = '\Users\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\Testing the new mask';

%Write name of the original mask here.
MaskFileName = 'tenframemaskFromPNG.tif';

%Write the name of the ND2 file to overlay here
BFFileName = 'ana33047_nh4_to_minusn_part1_0002_crop_crop(every10th)_shortened.nd2';

%Write what channel number to overlay for Fluor here
ChNum = 2;

%Write the desired name for the overlay mask (no numbers) here. 
%Must be the same name for output and input functions
maskName = 'maskOverlay10';

%If the mask starts some number of frames in, enter the difference from 1
%here
CountDiff = 0;

%Write the directory for the output here.
outDir = '\Users\Chris\Documents\Boulder\Lab Work\MatLab\Anabaena\Testing the new mask/Export';

%% Make the Movie %%

%First, load the metadata from the mask and start getting data from
%brightfield

MaskInfo = imfinfo(fullfile(currDir, MaskFileName));
BFFile =  BioformatsImage ( fullfile(currDir, BFFileName));

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

        % Overlay the mask as Blue and the brightfield as red and blue
        currOverlay = cat(3, currBF, currMask * .5, currBF);

