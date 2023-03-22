          % % --- Segment Fremyella --- % %
          

          
            % % --- User Inputs --- % %
            
            clear
            
            Directory = 'E:\Movies\first movie\';
            MovieName = '2020.3.5_ana33047_minusn_0003_Crop.nd2';
            
            ExportFile = 0; %1 to export a csv file with data, 0 to skip
            ExportName = 'Export_2020.12.8_fremyella_filament_to_hormogonia_0002.csv';
            
            PlayMovie = 0; %1 for yes, 0 for no
            FrameRate = 20; % This framterate does not include the time for your computer to run the script, and so the observed framerate will likely be slower
            
            ShowPlots = 0; %1 for yes, 0 for no
            
            % % --- User Inputs --- % %
            
            
            
            % % --- Segmentation and Data --- % %
            
            File =  BioformatsImage ( fullfile(Directory,MovieName));
            
            Data = struct([]);
            Results.AverageArea = [];
            Results.CellPixelOverBackground = [];
            Results.ChangeMass = [];
            Results.CellToCell = [];
            Results.BackgroundToCell = [];
            Results.CellToBackground = [];
            Results.NumObj = [];
            Results.ChangeObjNum = [];
            
            CurrCellPixels = [];
            LastCellPixels = [];
            
            BFLast = [];
            
            % for iFrame = 1:File.sizeT
            for iFrame = 70
                
                % % Segment down to single pixels % %
                currBF = getPlane(File,1,1,iFrame);
                currBFroi = currBF;
                currBFO = getPlane(File,1,2,iFrame);
                
                figure(1)
                clf('reset')
                imshow(currBFroi)
                roi = drawrectangle('StripeColor','y');
                rect = [roi.Vertices(1,1),roi.Vertices(1,2),(roi.Vertices(3,1)-roi.Vertices(1,1)),(roi.Vertices(2,2)-roi.Vertices(1,2))];
                
                currBFO = imcrop(currBFO,rect);
%                 currBFO = 
                
                BF =  ~imbinarize(currBFO);
                BF = imdilate(BF,strel ('disk',6));
                BF = imerode(BF, strel ('rectangle',[8 2]));
                BF = ~BF;
                BFi = bwmorph(BF,'shrink',Inf);
                BFi = imdilate(BFi,strel ('disk',2));
                BFi = bwmorph(BFi,'shrink',Inf);
                
                % % Orient the pixels to align with the x axis % %
                
                % find a line between 1 and final and flatten them, use an average of Y
                % points to make the line?
                
                celldata = regionprops(BFi,'Centroid');
                
                X = NaN (1, numel(celldata));
                Y = NaN (1, numel(celldata));
                
                for iCell = 1:numel(celldata)
                    
                    X(iCell) = celldata(iCell).Centroid(1);
                    Y(iCell) = celldata(iCell).Centroid(2);
                end
                
                coefficients = polyfit(X, Y, 1);
                yFitted = polyval(coefficients, X);
                
                hold off
                hold on
                %
                %      imshow(BF)
                %
                %      plot(X,yFitted)
                
                hold off
                
                slope = (yFitted(numel(celldata)) - yFitted(1)) ./ (X(numel(celldata)) - X(1));
                
                angle = atand(slope);
                
                BF = imrotate(BF,angle,'bilinear','loose');
                
                BF = bwmorph(BF,'shrink',Inf);
                BF = imdilate(BF,strel ('disk',2));
                BF = bwmorph(BF,'shrink',Inf);
                
                celldataR = regionprops(BF,'Centroid');
                
                X = NaN (1, numel(celldataR));
                Y = NaN (1, numel(celldataR));
                
                for iCell = 1:numel(celldataR)
                    
                    X(iCell) = celldataR(iCell).Centroid(1);
                    Y(iCell) = celldataR(iCell).Centroid(2);
                end
                
                coefficients = polyfit(X, Y, 1);
                yFitted = polyval(coefficients, X);
                meanY = mean (yFitted);
                
                for iCell = 1:numel(celldataR)
                    
                    celldataR(iCell).Centroid(2) = celldataR(iCell).Centroid(2) - meanY;
                    
                end
                
                for iCell = 1:numel(celldataR)
                    
                    Y(iCell) = -celldataR(iCell).Centroid(2);
                end
                
                yu = max(Y);
                yl = min(Y);
                yr = (yu-yl);                               % Range of ‘y’
                yz = Y-yu+(yr/2);
                zx = X(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
                per = 2*mean(diff(zx));                     % Estimate period
                ym = mean(Y);                               % Estimate offset
                
                fit = @(b,X)  (b(1).*(sin(2*pi*X./b(2) + 2*pi/b(3))) + b(4));    % Function to fit
                fcn = @(b) sum((fit(b,X) - Y).^2);                              % Least-Squares cost function
                s = fminsearch(fcn, [yr;  per;  -1;  ym]);                       % Minimise Least-Squares
                
                xp = linspace(min(X),max(X));
                
                imshow(BF)
                
                
                figure(2)
                clf('reset')
                plot(X,Y,'b',  xp,fit(s,xp), 'r')
                grid
                
                currBFr = imcrop(currBF,rect);
                currBFr = imrotate(currBFr,angle,'bilinear','loose');
                
                figure(3)
                clf('reset')
                hold on
                imshow(currBFr)
                plot( xp,-(fit(s,(xp))-meanY), 'r')
                plot(X,-(Y-meanY),'b')
                hold off
                
                
                
                
                
                





















%     celldata = regionprops(BF,'Area','Centroid','PixelList');
    
    
%     Data{iFrame} = celldata;
%     
%     Results.AverageArea(iFrame) = mean ([celldata.Area]);
%     Results.CellPixelOverBackground(iFrame) = sum([celldata.Area]) / ((File.width * File.height) - sum([celldata.Area]));
%     
%     
%     if iFrame == 1;
%         
%         BFLast = BF;
%         
%         Results.ChangeMass(iFrame) = 0;
%         Results.CellToCell(iFrame) = 0;
%         Results.BackgroundToCell(iFrame) = 0;
%         Results.CellToBackground(iFrame) = 0;
%         
%         Results.NumObj(iFrame) = numel(size(celldata.Area))+1; %why does calling the size give you the size minus 1?
%         Results.ChangeObjNum(iFrame) = 0;
%         
%     else
%         
%         CellMap = 2 * BF + BFLast;
%         Cells = CellMap == 3;
%         GCells = CellMap == 2;
%         LCells = CellMap == 1;
%         
%         Cellsdata = regionprops(Cells,'Area');
%         GCellsdata = regionprops(GCells,'Area');
%         LCellsdata = regionprops(LCells,'Area');
%         
%         Results.ChangeMass(iFrame) = Results.CellPixelOverBackground(iFrame)-Results.CellPixelOverBackground(iFrame-1);
%         Results.CellToCell(iFrame) = sum([Cellsdata.Area]);
%         Results.BackgroundToCell(iFrame) = sum([GCellsdata.Area]);
%         Results.CellToBackground(iFrame) = sum([LCellsdata.Area]);
%         
%         BFLast = BF;
%         
%         Results.NumObj(iFrame) = numel(size(celldata.Area))+1; %why does calling the size give you the size minus 1?
%         Results.ChangeObjNum(iFrame) = Results.NumObj(iFrame) - Results.NumObj(iFrame-1);
%         
%     end
%     
%     if PlayMovie == 1;
%         
%         figure(1);
%         imshowpair(bwperim(BF),currBF)
%         pause(1/FrameRate);
%         
%     elseif PlayMovie == 0;
%         
%     else
%         
%         error('Play Movie Invalid Value, Enter 1(Yes) or 0(No)')
%         
%     end
    
end


%             % % --- Plot Results --- % %
% 
% if ShowPlots == 1;
%     
%     figure(1);
%     plot(1:File.sizeT,[Results.AverageArea]);
%     title('Average Area');
%     
%     figure(2);
%     plot(1:File.sizeT,[Results.CellPixelOverBackground]);
%     title('Cell Area Over Total Area');
%     
%     figure(3)
%     plot(1:File.sizeT,[Results.ChangeMass]);
%     title('Change in Cell Mass');
%     
%     figure(4)
%     plot(1:File.sizeT,[Results.CellToCell]);
%     title('Cell To Cell');
%     
%     figure(5)
%     plot(1:File.sizeT,[Results.BackgroundToCell]);
%     title('Background To Cell');
%     
%     figure(6)
%     plot(1:File.sizeT,[Results.CellToBackground]);
%     title('Cell To Background');
%     
%     figure(7)
%     plot(1:File.sizeT,[Results.NumObj]);
%     title('Number of Objects');
%     
%     figure(8)
%     plot(1:File.sizeT,[Results.ChangeObjNum]);
%     title('Change in Number of Objects');
%     
% elseif ShowPlots == 0;
%     
% else
%     
%     error('Show Plots Invalid Value, Enter 1(Yes) or 0(No)')
%     
% end
% 
% 
%             % % --- Export csv --- % %
% 
% if ExportFile == 1;
%     
%     T = [];
%     
%     T(1,:) = 1:File.sizeT;
%     T(2,:) = [Results.AverageArea];
%     T(3,:) = [Results.CellPixelOverBackground];
%     T(4,:) = [Results.ChangeMass];
%     T(5,:) = [Results.CellToCell];
%     T(6,:) = [Results.BackgroundToCell];
%     T(7,:) = [Results.CellToBackground];
%     T(8,:) = [Results.NumObj];
%     T(9,:) = [Results.ChangeObjNum];
%     T = T';
%     
%     T = array2table(T,'VariableNames',{'Frame Number','Average Area',...
%         'Cell Pixels Over Background','Change in Cell Mass','Cell To Cell',...
%         'Background To Cell','Cell To Background','Number of Objects','Change in Number of Objects'});
%     
%     
%     writetable(T,fullfile(Directory,ExportName));
%     
% % elseif Export == 0;
%     
% else
%     
%     error('Export Invalid Value, Enter 1(Yes) or 0(No)')
%     
% end
% 
% 
% clearvars 
% clc
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
