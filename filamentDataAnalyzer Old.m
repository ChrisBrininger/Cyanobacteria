classdef filamentDataAnalyzer < DataAnalyzer
    %DATAANALYZER  Data analysis for cyanobacteria cells
    %
    %  OBJ = DATAANALYZER creates an empty DataAnalyzer object. Use the
    %  method importdata to load TrackArray data into the object.
    
    properties
        
    end
    
    methods
        
        function obj = filamentDataAnalyzer(varargin)
            %FILAMENTDATAANALYZER  Construct a new object
            %
            %  OBJ = FILAMENTDATAANALYZER will create an empty object. To
            %  import data use IMPORTDATA(OBJ).
            %
            %  OBJ = FILAMENTDATAANALYZER(F) will call importData on the
            %  file specified. F should be a string containing full path
            %  names to MAT-file for import.
            %
            %  See also: IMPORTDATA
            
            if ~isempty(varargin)
                obj = importdata(obj, varargin{:});
            end
            
        end
        
        function obj = analyze(obj)
            
            obj = analyze@DataAnalyzer(obj);
            
            %~~~ Then, find the filament positions ~~~%
            %NOTE: This will currently break if the cellID of
            %filamentPosition 1 divides. It will also break if any of the
            %cell's divide on the second frame.
            
            %First, initialize the field filamentPositions in the data
            %struct for each cellID, of the correct length.
            for iCell = 1:numel(obj.Tracks)
                
                numFrames = numel(obj.Tracks(iCell).Frames);
                obj.Tracks(iCell).Data.FilamentPosition = cell(1, numFrames);
                
            end
            
            %Then, find and fill filamentPosition by frame. For first
            %frame, use getLineData function to initialize filamentPosition
            %for each cell. For all subsequent frames, fill in each cell's
            %filament position, starting from 1, based on its filament
            %position from the previous frame, and from previous cells of
            %the current frame, for division events.
            for iFrame = 1:numel(obj.FileMetadata.Timestamps)
                
                if iFrame == 1 %use Centroid positions to find filamentPosition
                    
                    centroidCoords = [];
                    currcellIDs = [];
                    for iCell = 1:numel(obj.Tracks)
                        
                        if obj.Tracks(iCell).Generation == 1
                            centroidCoords(end+1, 1) = obj.Tracks(iCell).Data.Centroid{1}(1);
                            centroidCoords(end, 2) = obj.Tracks(iCell).Data.Centroid{1}(2);
                            currcellIDs(end+1) = obj.Tracks(iCell).ID;
                        end
                        
                    end
                    roundedCoords = round(centroidCoords);
                    
                    S = filamentDataAnalyzer.getLineData(roundedCoords);
                    
                    currfilPosits = [];
                    for iCell = currcellIDs
                        
                        currXCoord = roundedCoords(iCell, 1);
                        currYCoord = roundedCoords(iCell, 2);
                        indX = find(ismember(S.SortedCoords(:, 1), currXCoord));
                        indY = find(ismember(S.SortedCoords(:, 2), currYCoord));
                        %In case multiple centroids share the same X or Y
                        %coordinate, find the shared index
                        commonInd = intersect(indX, indY);
                        
                        obj.Tracks(iCell).Data.FilamentPosition{1} = commonInd;
                        
                        %For easily keeping track of filament position of
                        %known cells
                        currfilPosits(end+1) = commonInd;
                        
                    end
                    
                else %If iFrame ~= 1, find filamentPosition based off the previous frame
                    
                    %First, find the cellID from previous frame of
                    %filamentPosition 1.
                    nextCellIDs = []; %fill throughout this frame, for preparation for next frame.
                    nextfilPosits = []; %fill throughout this frame, for preparation for next frame.
                    shiftNumNextFrame = 0; %Counter if multiple cells divide in same frame.
                    skipCell = 0; %Used to skip analysis of 2nd daughter cell, as positions of both are determined simultaneously.
                    for iCell = 1:numel(currfilPosits)
                        
                        if skipCell == 1
                            %Skip analysis, as was performed in previous
                            %iCell, due to new division event.
                            skipCell = 0;
                            
                        else
                            
                            currFilPos = currfilPosits(iCell);
                            currID = currcellIDs(iCell);
                            newDivision = 0; %Required for properly giving ID and position info in preparation for next frame.
                            
                            %This if statement asks if we are dealing with
                            %a newly divided cell or not. If nan, then it
                            %is a newly divided cell.
                            if ~isnan(currFilPos)
                                %filamentPosition is same as in last frame.
                                obj.Tracks(currID).Data.FilamentPosition{find(obj.Tracks(currID).Frames == iFrame)} = currFilPos;
                                
                            else
                                %Assume there is another nan right after
                                %this one (for the other daughter).
                                %Determine each daughters centroid position
                                %relative to the cellID of the previous
                                %cell (furthest determined in the
                                %filament).
                                
                                newDivision = 1;
                                
                                d1Centroid = obj.Tracks(currID).Data.Centroid{1};
                                
                                otherDaughterID = currcellIDs(iCell+1);
                                d2Centroid = obj.Tracks(otherDaughterID).Data.Centroid{1};
                                
                                prevCellID = currcellIDs(iCell-1);
                                prevCellCentroid = obj.Tracks(prevCellID).Data.Centroid{find(obj.Tracks(prevCellID).Frames == iFrame)};
                                
                                %Calculate distances
                                dist1 = sqrt((d1Centroid(1)-prevCellCentroid(1)).^2+(d1Centroid(2)-prevCellCentroid(2)).^2);
                                dist2 = sqrt((d2Centroid(1)-prevCellCentroid(1)).^2+(d2Centroid(2)-prevCellCentroid(2)).^2);
                                
                                if dist1 < dist2
                                    %In this case, cellIDs are in correct
                                    %order.
                                    
                                    %Replace nans (required if two cells next
                                    %to each other divide on same frame.)
                                    currfilPosits(iCell) = max(currfilPosits(1:iCell))+1;
                                    currfilPosits(iCell+1) = currfilPosits(iCell)+1;
                                    
                                    obj.Tracks(currID).Data.FilamentPosition{find(obj.Tracks(currID).Frames == iFrame)} = currfilPosits(iCell);
                                    obj.Tracks(otherDaughterID).Data.FilamentPosition{find(obj.Tracks(currID).Frames == iFrame)} = currfilPosits(iCell)+1;
                                    
                                else
                                    %cellIDs are not in correct order.
                                    %Switch position of cellIDs so that the
                                    %filament order is still in ascending
                                    %order.
                                    currcellIDs(iCell) = currcellIDs(iCell+1);
                                    currcellIDs(iCell+1) = currID;
                                    currID = currcellIDs(iCell);
                                    
                                    %Then, replace nans.
                                    currfilPosits(iCell) = max(currfilPosits(1:iCell))+1;
                                    currfilPosits(iCell+1) = currfilPosits(iCell)+1;
                                    
                                    obj.Tracks(currID).Data.FilamentPosition{find(obj.Tracks(currID).Frames == iFrame)} = currfilPosits(iCell);
                                    obj.Tracks(currcellIDs(iCell+1)).Data.FilamentPosition{find(obj.Tracks(currID).Frames == iFrame)} = currfilPosits(iCell+1);
                                    
                                end
                                
                            end
                            
                            %This if statement checks if we are in the
                            %final frame of the current ID. If so, we
                            %replace this ID for the next frame with its
                            %two daughters, with unknown (nan)
                            %filamentPositions. This is in preparation for
                            %the next iFrame iteration, and does not impact
                            %the current frame.
                            if obj.Tracks(currID).Frames(find(obj.Tracks(currID).Frames == iFrame)) == obj.Tracks(currID).Frames(end)
                                
                                %If the cell is going to divide next frame,
                                %replace its ID in cellIDs with its two
                                %daughter cells, with unknown filamentPositions
                                nextCellIDs = [nextCellIDs, obj.Tracks(currID).DaughterID];
                                nextfilPosits = [nextfilPosits, nan, nan];
                                
                                %If cell divides, all future positions must
                                %increase by one due to presence of an
                                %additional cell. This corrects for
                                %frameshift required if multiple cells
                                %divide in same frame.
                                shiftNumNextFrame = shiftNumNextFrame + 1;
                                
                            else
                                %If current cell is not in the final
                                %frame, must check if current cell is in
                                %its first frame (newDivision). If so, the
                                %filamentPositions of both daughters were
                                %determined simultaneously, so both
                                %daughters IDs and filPositions must be
                                %incorporated into the next frame's IDs and
                                %positions.
                                if newDivision == 0
                                    %Cell is not newly divided (in its
                                    %first frame of existence). Simply
                                    %append its Id and position for
                                    %preparation for next frame.
                                    nextCellIDs = [nextCellIDs, currID];
                                    nextfilPosits = [nextfilPosits, currfilPosits(iCell) + shiftNumNextFrame];
                                    
                                    skipCell = 0;
                                else
                                    %cell is newly divided. Add info from
                                    %both cells for preparation for next
                                    %frame.
                                    nextCellIDs = [nextCellIDs, currcellIDs(iCell), currcellIDs(iCell+1)];
                                    nextfilPosits = [nextfilPosits, currfilPosits(iCell) + shiftNumNextFrame, currfilPosits(iCell+1) + shiftNumNextFrame];
                                    
                                    %Then, Skip the next icell, since both daughters
                                    %determined at once
                                    skipCell = 1;
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                    %Check that we are not in final frame.
                    if iFrame ~= numel(obj.FileMetadata.Timestamps)
                        %Reset the 'next' IDs to the current IDs, in
                        %preparation for the next frame
                        currcellIDs = nextCellIDs;
                        currfilPosits = nextfilPosits;
                    end
                    
                end
                
            end
            %~~~~End of finding FilamentPosition~~~~%
            
            
            
            % % --- Fix growth rates --- % %
            
            for ii = 1:numel(obj.Tracks)
                
                %--- Calculate Growth Rate ---%
                tt = (obj.Tracks(ii).Frames) * obj.MeanDeltaT;
                
                %Replace empty values (for skipped frames) with NaNs
                len = (([obj.Tracks(ii).Data.MajorAxisLength{:}]+[obj.Tracks(ii).Data.MinorAxisLength{:}])-16);
                
                %                 idxEmpty = find(cellfun(@isempty, len));
                %                 for iC = 1:numel(idxEmpty)
                %                     len{idxEmpty(iC)} = NaN;
                %                 end
                %
                len = len * obj.FileMetadata.PhysicalPxSize(1);
                
                [obj.Tracks(ii).GrowthRate, obj.Tracks(ii).GRFitY, obj.Tracks(ii).GRFitRes] = ...
                    DataAnalyzer.fitGrowthRate(tt, len);
                
                %Put fields in useful order, followed by everything else in
                %alphabetical order
                order = {'ID', 'Frames', 'Data', 'MotherID', 'DaughterID', 'GrowthRate',...
                    'Generation', 'GRFitY', 'GRFitRes'};
                order = order';
                
                fieldOrder = fieldnames(obj.Tracks);
                fieldOrder(ismember(fieldOrder, order)) = [];
                fieldOrder = sort(fieldOrder);
                fieldOrder = [order; fieldOrder];
                
                obj.Tracks = orderfields(obj.Tracks, fieldOrder);
                
            end
            
            obj = findheterocysts(obj);
            
        end
        
        function obj = findheterocysts (obj)
            
            for ii = 1:numel(obj.Tracks)
                
                %--- Calculate Growth Rate ---%
                tt = (obj.Tracks(ii).Frames) * obj.MeanDeltaT;
                
                %Replace empty values (for skipped frames) with NaNs
                len = (([obj.Tracks(ii).Data.MajorAxisLength{:}]+[obj.Tracks(ii).Data.MinorAxisLength{:}])-16);
                
                %                 idxEmpty = find(cellfun(@isempty, len));
                %                 for iC = 1:numel(idxEmpty)
                %                     len{idxEmpty(iC)} = NaN;
                %                 end
                %
                len = len * obj.FileMetadata.PhysicalPxSize(1);
                
                [obj.Tracks(ii).GrowthRate, obj.Tracks(ii).GRFitY, obj.Tracks(ii).GRFitRes] = ...
                    filamentDataAnalyzer.fitGrowthRate(tt, len);
                
                %--- Find Heterocysts ---%
                
                endInt = ([obj.Tracks(ii).Data.TotalIntCy5{end}]/numel(obj.Tracks(ii).Data.RegisteredPxInd{end}))...
                    /([obj.Tracks(ii).Data.TotalIntCy5{1}]/numel(obj.Tracks(ii).Data.RegisteredPxInd{1}));
                
                if endInt <=.80
                    
                    obj.Tracks(ii).isHeterocyst=1;
                    
                else
                    
                    obj.Tracks(ii).isHeterocyst=0;
                    
                end
            end
            
        end
        
        
        
        
        function obj = NearestHeterocyst (obj)
            
            % % set all heterocyst identities to 0 to clear the slate% %
            for ii = 1:numel(obj.Tracks)
                obj.Tracks(ii).isHeterocyst=0;
            end
            
            % % Start by analyzing the movie one frame at a time. % %
            for iFrame = 1:numel(obj.FileMetadata.Timestamps)
                
                ListofHet=[];
                
                HetDirection=nan;
                
                % see if any cells are heterocysts in this frame
                for ii = 1:numel(obj.Tracks)
                    
                    % does this cell exist this frame? if not, skip
                    
                    E = find(obj.Tracks(ii).Frames == iFrame);
                    
                    if isempty(E) == [0]
                        
                        % if the cell is not currently a heterocyst, check to see
                        % if it should be
                        if obj.Tracks (ii).isHeterocyst == 0
                            
                            
                            
                            %--- Find Heterocysts ---%
                            
                            CurrInt = ([obj.Tracks(ii).Data.TotalIntCy5{E}]/numel(obj.Tracks(ii).Data.RegisteredPxInd{E}))...
                                /([obj.Tracks(ii).Data.TotalIntCy5{1}]/numel(obj.Tracks(ii).Data.RegisteredPxInd{1}));
                            
                            %                          Current Heterocyst threshold set at .8 of
                            %                          initial intensity
                            if CurrInt <=.80
                                
                                obj.Tracks(ii).isHeterocyst=1;
                                
                                
                            end
                            
                        end
                        
                        %assign heterocyst cell numbers to a list
                        if obj.Tracks (ii).isHeterocyst == 1
                            
                            ListofHet(end+1) = ii;
                            
                        end
                    end
                end
                
                
                
                % Determine the distance of each non heterocyst to its
                % closest heterocyst
                if isempty(ListofHet) == 0
                    
                    % for each cell again
                    for ii = 1:numel(obj.Tracks)
                        
                        Distance = [];
                        
                        
                        % Find out if it exist during the current frame
                        E = find(obj.Tracks(ii).Frames == iFrame);
                        
                        if isempty(E) == [0]
                            
                            if obj.Tracks(ii).isHeterocyst == 0
                                
                                for Hets = (ListofHet)
                                    
                                    %find out the correct cell frame for
                                    %the heterocyst
                                    Eh = find(obj.Tracks(Hets).Frames == iFrame);
                                    
                                    Distance(end+1) = obj.Tracks(Hets).Data.FilamentPosition{Eh}...
                                        - obj.Tracks(ii).Data.FilamentPosition{E};
                                    
                                end
                                
                                % take the absolute value, then the minimum
                                % value to find the smallest cell count
                                % difference and report it
                                aDistance = abs(Distance);
                                
                                mDistance = min(aDistance);
                                if ii == 79
                                    if E==22
                                        keyboard
                                    end
                                end
                                obj.Tracks(ii).Data.CellsToHeterocysts (E) = mDistance;
                                
                                
                                
                                % % --------------------------------------------------------------------% %
                                % This starts the code looking to find the uM
                                DistCellNum = [];
                                
                                % Create a blank list for the cell
                                % IDs
                                DistCells=[];
                                DistCellsP=[];
                                DistCellsN=[];
                                
                                if isempty(find(Distance==mDistance))==1
                                    HetDirection = -1;
                                    % Finds the distance in uM if the
                                    % CellToHeterocysts is Negative
                                    %                                        keyboard
                                else
                                    if numel(find(mDistance==aDistance))==1
                                        HetDirection = 1;
                                        % Finds the distance in uM if the
                                        % CellToHeterocysts is Positive
                                        % and there is only 1 value
                                        
                                    else
                                        HetDirection = 2;
                                        %Reserving this space for a case
                                        %in which there are 2
                                        %heterocysts, one + one -, that
                                        %are the same count of cells
                                        %away
                                        %                                            error('HetDirection = 2')
                                    end
                                end
                                
                                if HetDirection <= 1
                                    % Finds what the filament positions
                                    % to the direction of the cell are
                                    % for the number cells between it
                                    % and the heterocyst (runs both
                                    % directions for + and -)
                                    DistCellFillNum = [(obj.Tracks(ii).Data.FilamentPosition{E}+HetDirection):(obj.Tracks(ii).Data.FilamentPosition{E}+((mDistance-1)*HetDirection)),...
                                        (obj.Tracks(ii).Data.FilamentPosition{E}+((mDistance-1)*HetDirection)):(obj.Tracks(ii).Data.FilamentPosition{E}+HetDirection)];
                                    
                                    % check each cell to see if it is in
                                    % one of the filament positions in
                                    % this frame
                                    for eCell = 1:numel(obj.Tracks)
                                        
                                        % Check to see if cell exists in
                                        % this frame
                                        Ec = find(obj.Tracks(eCell).Frames == iFrame);
                                        
                                        DistCell = [];
                                        
                                        if isempty(Ec)==0;
                                            
                                            % If it does, check to see if
                                            % it has a correct filament
                                            % position in that frame
                                            DistCell = find (obj.Tracks(eCell).Data.FilamentPosition {Ec}==DistCellFillNum);
                                            
                                            if isempty(DistCell)==0
                                                %                                             % if a cell was found,
                                                %                                             add its length to
                                                %                                             DistCells list
                                                % Cell Length Calculation
                                                CellLength = obj.Tracks(eCell).Data.MajorAxisLength {Ec}+obj.Tracks(eCell).Data.MinorAxisLength {Ec}-16;
                                                
                                                DistCells(end+1) = CellLength;
                                                
                                            end
                                        end
                                    end
                                else
                                    %This does the same thing as above
                                    %twice, assuming a + and -
                                    %direction, and then uses the
                                    %smaller value to report
                                    %                                        keyboard
                                    % Finds what the filament positions
                                    % to the direction of the cell are
                                    % for the number cells between it
                                    % and the heterocyst
                                    DistCellFillNumP = (obj.Tracks(ii).Data.FilamentPosition{E}+1):(obj.Tracks(ii).Data.FilamentPosition{E}+mDistance-1);
                                    DistCellFillNumN = (obj.Tracks(ii).Data.FilamentPosition{E}-mDistance+1):(obj.Tracks(ii).Data.FilamentPosition{E}-1);
                                    
                                    % check each cell to see if it is in
                                    % one of the filament positions in
                                    % this frame
                                    for eCell = 1:numel(obj.Tracks)
                                        
                                        % Check to see if cell exists in
                                        % this frame
                                        Ec = find(obj.Tracks(eCell).Frames == iFrame);
                                        
                                        DistCellP = [];
                                        DistCellN = [];
                                        
                                        if isempty(Ec)==0;
                                            
                                            % If it does, check to see if
                                            % it has a correct filament
                                            % position in that frame
                                            DistCellP = find (obj.Tracks(eCell).Data.FilamentPosition {Ec}==DistCellFillNumP);
                                            DistCellN = find (obj.Tracks(eCell).Data.FilamentPosition {Ec}==DistCellFillNumN);
                                            
                                            if isempty(DistCellP)==0
                                                %                                                 % if a cell was found,
                                                %                                                 add its length to
                                                %                                                 DistCells list
                                                CellLength = obj.Tracks(eCell).Data.MajorAxisLength {Ec}+obj.Tracks(eCell).Data.MinorAxisLength {Ec}-16;
                                                
                                                DistCellsP(end+1) = CellLength;
                                            else
                                                if isempty(DistCellN)==0
                                                    CellLength = obj.Tracks(eCell).Data.MajorAxisLength {Ec}+obj.Tracks(eCell).Data.MinorAxisLength {Ec}-16;
                                                    
                                                    DistCellsN(end+1) = CellLength;
                                                end
                                            end
                                        end
                                    end
                                    
                                    DistCellsP = sum(DistCellsP);
                                    DistCellsN = sum(DistCellsN);
                                    Distances = [DistCellsP,DistCellsN];
                                    DistCells = min(Distances);
                                    
                                end
                                %Sum up all dist cells found and
                                %make that value the  microns
                                obj.Tracks(ii).Data.MicronsToHeterocysts (E) = sum(DistCells);
                                
                            end
                        end
                    end
                end
                
            end
            
            % Make it so that cells that never coexist with a heterocyst still
            % have the fields and fill them with 0s for graphing
            for iCell =  1:numel(obj.Tracks)
                
                if isfield (obj.Tracks(iCell).Data,'CellsToHeterocysts') == 0
                    
                    for eFrame = 1:numel(obj.Tracks(iCell).Frames)
                        
                        obj.Tracks(iCell).Data.CellsToHeterocysts (eFrame) = 0;
                        
                    end
                end
                
                if isfield (obj.Tracks(iCell).Data,'MicronsToHeterocysts') == 0
                    
                    for eFrame = 1:numel(obj.Tracks(iCell).Frames)
                        
                        obj.Tracks(iCell).Data.MicronsToHeterocysts (eFrame) = 0;
                        
                    end
                end
                
            end
        end
        
        
        function obj = fixTimeStamps(obj, deltaTinHours)
            %FIXTIMESTAMPS manually fixes time stamps of nd2 files. Only
            %required when the nd2 file timestamps are clearly incorrect
            %(extremely large or small numbers, or zero). This is sometimes a
            %problem when using a cropped (xy) nd2 file.
            %
            %  FIXTIMESTAMPS(obj, deltaTinHours) fixes the timestamps given
            %  the deltaTinHours. The first TimeStamp is set to zero hours.
            %  All values in tblData (ie GrowthRate) are then recalculated
            %  using the corrected TimeStamp data.
            
            for iFID = 1:numel(obj.FileMetadata)
                obj.FileMetadata.Timestamps = (1:numel(obj.FileMetadata.Timestamps)) .* deltaTinHours;
                obj.FileMetadata.Timestamps = (obj.FileMetadata.Timestamps - obj.FileMetadata.Timestamps(1));% .* 3600;
                obj.FileMetadata.MeanDeltaT = deltaTinHours;
                obj.MeanDeltaT = deltaTinHours;
            end
            
            %re-analyze the data
            obj = analyze(obj);
            
        end
        
        
        %~~~~~~  Plotting Functions  ~~~~~~%
        
        function plotLineage(obj, cellID, varargin)
            %PLOTLINEAGE  Plot cell lineage
            %Set to add axes and subtract a previously measure width
            %constant
            %  PLOTLINEAGE(OBJ, ID) plots the cell length over time for
            %  cell ID and all its ancestors.
            
            %Parse the variable input argument
            plotfit = false;
            while ~isempty(varargin)
                if ischar(varargin{1}) && strcmpi(varargin{1}, 'plotfit')
                    %Plot fitted growth rate
                    plotfit = true;
                end
                varargin(1) = [];
            end
            
            while ~isnan(cellID)
                
                %                 Randomly Selected Cell Widths (px)
                %                   16, 15, 16.7,
                
                tVec = obj.FileMetadata.Timestamps(obj.Tracks(cellID).Frames);
                lenInMicrons = (([obj.Tracks(cellID).Data.MajorAxisLength{:}]+[obj.Tracks(cellID).Data.MinorAxisLength{:}])-16) .* obj.FileMetadata.PhysicalPxSize(1);
                
                lineWidth = 2.5;
                markerSize = 1;
                
                plot(tVec, lenInMicrons, 'b.-', 'MarkerSize', markerSize, 'LineWidth', lineWidth)
                hold on
                
                if plotfit
                    m = obj.Tracks(cellID).GrowthRate;
                    yint = obj.Tracks(cellID).GRFitY;
                    plot(tVec, exp(m * tVec + yint),'k-', 'LineWidth', lineWidth-1)
                    
                end
                
                cellID = obj.Tracks(cellID).MotherID;
                
            end
            
            %             Graph Aesthetics
            
            hold off
            %             xlim([0 65])
            xlabel('Time (hours)')
            ylabel('Cell length (\mum)')
            set(findall(gcf,'-property','FontSize'),'FontSize',17)
            set(gca,'linewidth',0.5)
            box off
        end
    end
    
    methods (Static)
        
        function lineData = getLineData(coords)
            %GETLINEDATA  Measures data of a line
            %
            %  S = GETLINEDATA(C) returns a struct S containing data about
            %  the line specified by the coordinate vector C. C must be an
            %  Nx2 list of coordinates where C(:, 1) is x and C(:, 2) is y.
            %  It is expected that the points in C is connected by no more
            %  than 1.4 pixels (i.e. the object mask is 8-connected).
            %
            %  S has the following fields:
            %     Centroid = Center coordinate of the line
            %     Length = Length of the line
            %     SortedCoords = Sorted coordinates of the line
            %     Excluded = Any coordinates that were excluded
            %
            %  The SortedCoords will be an Nx2 matrix specifying
            %  coordinates going from one end of the line to the other.
            %
            %  Excluded contains a list of coordinates that were excluded
            %  from the line. These are points that did not connect to the
            %  traced line, e.g. if there was a branch.
            %
            %  The algorithm works by first tracing the points to find an
            %  end point. The cumulative distance to the end point is then
            %  computed and used to order the points as well as getting the
            %  length of the line. The middle point is then the point that
            %  is closest to the length/2.
            %
            %  Example:
            %  %Assume after segmentation and skeletonization you have a
            %  %mask. Note the use of 'PixelList' instead of
            %  %'PixelIdxList'.
            %  data = regionprops(mask, 'PixelList');
            %
            %  lineData = ActinTracker.findLineCenter(data.PixelList);
            
            
            %Initialize variables
            isSorted = [true; false(size(coords, 1) - 1, 1)];  %Flag if point has been sorted
            sortIdx = [1; zeros(size(coords, 1) - 1, 1)];  %To store sorted indices
            ptrLastIndex = 1;  %Position of sort index to add to
            
            %Find an end point by travelling in one direction from the
            %first pixel.
            minInd = 1;  %Start at the first coordinate
            while ~all(isSorted)
                
                %Find the next nearest pixel
                sqDistToPtr = sum((coords - coords(minInd, :)).^2, 2);
                sqDistToPtr(isSorted) = Inf;
                
                [minDist, minInd] = min(sqDistToPtr);
                
                if minDist <= 50000
                    isSorted(minInd) = true;
                    
                    %Append to sorted indices
                    ptrLastIndex = ptrLastIndex + 1;
                    sortIdx(ptrLastIndex) = minInd;
                    
                else
                    break;
                end
            end
            
            %Shift the indices to the end of the array
            sortIdx = circshift(sortIdx, nnz(~isSorted));
            
            %Pointer will now count upwards so update the value
            ptrLastIndex = nnz(~isSorted) + 1;
            
            minInd = 1; %Reset the point back to the first coordinate
            while ~all(isSorted)
                
                %Find the next nearest pixel
                sqDistToPtr = sum((coords - coords(minInd, :)).^2, 2);
                sqDistToPtr(isSorted) = Inf;
                
                [minDist, minInd] = min(sqDistToPtr);
                
                if minDist <= 50000
                    isSorted(minInd) = true;
                    
                    %Add to sorted indices going upwards
                    ptrLastIndex = ptrLastIndex - 1;
                    sortIdx(ptrLastIndex) = minInd;
                    
                else
                    break;
                end
                
            end
            
            %Sort the array
            lineData.SortedCoords = coords(sortIdx, :);
            
            %Compute the line length
            distFromEnd = [0; cumsum(sqrt(sum((diff(lineData.SortedCoords)).^2, 2)))];
            lineData.Length = distFromEnd(end);
            
            %Find the center coordinate of the line
            [~, midPtLoc] = min(abs(distFromEnd - lineData.Length/2));
            lineData.Centroid = lineData.SortedCoords(midPtLoc, :);
            
            %Report any excluded data points
            lineData.Excluded = coords(~isSorted, :);
            
        end
        
        
    end
    
end












