% % Input Commands % %

% % % Note on limitations, currently the tif file can max out in size, so
% the script will limit the tif video to 350 frames per movie.

%Write the folder directory here.
currDir = 'E:\Movies\20210418_Ana_-N_-Buffer_+CaCO3';

%Write name of the first movie here.
MovieOneFileName = 'channelbf,cy5,cfp,rfp_seq0000_0000.nd2';

%Write the name of the second movie here.
MovieTwoFileName = 'channelbf,cy5,cfp,rfp_seq0000_0000.nd2';

%Write the name of the third movie here.
MovieThreeFileName = 'channelbf,cy5,cfp,rfp_seq0000_0000.nd2';

%Write the name of the fourth movie here.
MovieFourFileName = 'channelbf,cy5,cfp,rfp_seq0000_0000.nd2';

%Write the name of the fifth movie here.
MovieFiveFileName = 'channelbf,cy5,cfp,rfp_seq0000_0000.nd2';

%Write the name of the sixth movie here.
MovieSixFileName = 'channelbf,cy5,cfp,rfp_seq0000_0000.nd2';

%Write the desired name for the final movie (no numbers) here. 
%Must be the same name for output and input functions
movieName = 'CombinedMovie';

%Write the directory for the output here.
outDir = 'E:\Movies\20210418_Ana_-N_-Buffer_+CaCO3\CombinedMovieTesting\';

%Wirte the number of movies here. (up to 6 currently)
NumMovies = 6;

%Write the number of channels here.
NumChannels = 2;

%Write Scaling for Channel 1 Here (higher is brighter)
CH1Scale = 1;

%Write Scaling for Channel 2 Here (higher is brighter)
CH2Scale = 1;

%Write Scaling for Channel 3 Here (higher is brighter)
CH3Scale = 1;

%Write Scaling for Channel 4 Here (higher is brighter)
CH4Scale = 1;

%Write Scaling for Channel 5 Here (higher is brighter)
CH5Scale = 1;

%Write Scaling for Channel 6 Here (higher is brighter)
CH6Scale = 1;
%% 

Movie1Info = BioformatsImage(fullfile(currDir, MovieOneFileName));

Movie2Info = BioformatsImage(fullfile(currDir, MovieTwoFileName));

Movie3Info = BioformatsImage(fullfile(currDir, MovieThreeFileName));

Movie4Info = BioformatsImage(fullfile(currDir, MovieFourFileName));

Movie5Info = BioformatsImage(fullfile(currDir, MovieFiveFileName));

Movie6Info = BioformatsImage(fullfile(currDir, MovieSixFileName));

FrameCount1 = Movie1Info.sizeT;

FrameCount2 = Movie2Info.sizeT;

FrameCount3 = Movie3Info.sizeT;

FrameCount4 = Movie4Info.sizeT;

FrameCount5 = Movie5Info.sizeT;

FrameCount6 = Movie6Info.sizeT;

%Test
% 
% FrameCount1 = 9;
% 
% FrameCount2 = 9;
% 
% FrameCount3 = 9;
% 
% FrameCount4 = 9;
% 
% FrameCount5 = 9;
% 
% FrameCount6 = 9;

for iChannel = 1:NumChannels
    
    Channel = char(Movie1Info.channelNames(iChannel));
    
    if iChannel == 1;
        
        Scaling = CH1Scale;
        
    else if iChannel == 2;
            
            Scaling = CH2Scale;
            
        else if iChannel == 3;
                
                Scaling = CH3Scale;
                
            else if iChannel == 4;
                    
                    Scaling = CH4Scale;
                    
                else if iChannel == 5;
                        
                        Scaling = CH5Scale;
                        
                    else if iChannel == 6;
                            
                            Scaling = CH6Scale;
                            
                        else if iChannel >= 7;
                                
                                error('Too Many Channels')
                                
                            end
                        end
                    end
                end
            end
        end
    end
    
    if NumMovies >=7;
        
        error('Number of Movies Exceeds Current Limit (6)')
    end
    
    for iFrame = 1:FrameCount1;
        
        currFrame = getPlane(Movie1Info,1,iChannel,iFrame);
        
        if iFrame == 1
            imwrite((currFrame*Scaling), fullfile(outDir,[movieName,Channel,'.tif']), 'compression', 'none');
        else
            imwrite((currFrame*Scaling), fullfile(outDir,[movieName,Channel,'.tif']), 'writeMode', 'append', 'compression', 'none');
        end
    end
    
    if NumMovies >=2;
        
        for iFrame = 1:FrameCount2;
            
            currFrame = getPlane(Movie1Info,1,iChannel,iFrame);
            
            imwrite((currFrame*Scaling), fullfile(outDir,[movieName,Channel,'.tif']), 'writeMode', 'append', 'compression', 'none');
        end
        
        if NumMovies >=3;
            
            for iFrame = 1:FrameCount3;
                
                currFrame = getPlane(Movie1Info,1,iChannel,iFrame);
                
                imwrite((currFrame*Scaling), fullfile(outDir,[movieName,Channel,'.tif']), 'writeMode', 'append', 'compression', 'none');
            end
            
        end
        
        if NumMovies >=4;
            
            for iFrame = 1:FrameCount4;
                
                currFrame = getPlane(Movie1Info,1,iChannel,iFrame);
                
                imwrite((currFrame*Scaling), fullfile(outDir,[movieName,Channel,'.tif']), 'writeMode', 'append', 'compression', 'none');
            end
            
        end
        
        if NumMovies >=5;
            
            for iFrame = 1:FrameCount5;
                
                currFrame = getPlane(Movie1Info,1,iChannel,iFrame);
                
                imwrite((currFrame*Scaling), fullfile(outDir,[movieName,Channel,'.tif']), 'writeMode', 'append', 'compression', 'none');
            end
            
        end
        
        if NumMovies >=6;
            
            for iFrame = 1:FrameCount6;
                
                currFrame = getPlane(Movie1Info,1,iChannel,iFrame);
                
                imwrite((currFrame*Scaling), fullfile(outDir,[movieName,Channel,'.tif']), 'writeMode', 'append', 'compression', 'none');
            end
            
        end
        
    end
    
end