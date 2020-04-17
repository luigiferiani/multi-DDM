%{
% Version 1.3
% © Luigi Feriani, Maya Juenet, Pietro Cicuta, 2018 (luigi.feriani@gmail.com) 
% 
% DDM_Analysis.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% 
% Original work:
%
% Feriani, L., et al., Biohpysical Journal 2017
% "Assessing the Collective Dynamics of Motile Cilia in Cultures of Human
% Airway Cells"
%
% Chioccioli, M.*, Feriani, L.*, Kotar, J., Bratcher, P. E.**, Cicuta, P.**, Nature Communications 2019
% "Phenotyping ciliary dynamics and coordination in response to CFTR-modulators
% in Cystic Fibrosis respiratory epithelial cells"
% * These authors have contributed in equal measure
% ** These authors have contributed in equal measure
%
%}

classdef DDM_Analysis < matlab.mixin.Copyable
    %DDM_analysis class to handle and run ddm analysis on microscope videos
    %   Version 1.1.3 Changelog - 04/05/2018:
    %       changed way multiDDM is done, massive speedup if user has a
    %       gpu, good speedup without. 
    %       Methods affected:Variable_BoxSize_Analysis (but no change in
    %       interface with the user), DDM_on_Boxes became
    %       Setup_DDM_on_Boxes, added fit_Boxes
    %   Version 1.1.2 Changelog: 
    %       brand new motion detection algorithm, new properties
    %   Version 1.1.1 Changelog: 
    %       added support to 2 slightly different ways to find good boxes
    %       based on std_fs: by setting flag_ditch_imadjust it should work
    %       better on FOVs with just a few beating cells. This can be done
    %       a posteriori by calling gather_results with
    %       flag_recalculate_good_boxes set to true
    
    properties
        Filename            %Name of the video file
        Filepath            %Name of file path
        Height              %rows of video, (px)
        Width               %cols of video, (px)
        TotalNumberOfFrames %N_frames as written in the movie header
        FirstFrame          %first frame to analyse
        LastFrame           %last frame to analyse
        NumberOfFrames      %number of frames to analyse (LastFrame-FirstFrame+1)
        FrameRate           %
        Magnification       %for px2mum conversions
        std_fs              %standard deviation of pixel intensity through time
        lstd_sfd            % log10 of std of stdfilt differences between frames 2 apart
        detected_motion     % forreground objects from lstd_sfd
        Results             %results structure array
        %         fs                  %video as a 3d matrix. it will not be saved
        IsMovieCorrupted    %control variable for failsafe movie reading
        Kymograph
        Temperature         %temperature read from file
        SAVAlike            %structure with sava-like cbf measurements
    end
    
    properties (Hidden = true)
        px2mum              %pixels to micron conversion
        N_couple_frames_to_average %how many couples of frames to average in DDM_core
        max_tau
    end
    
    properties (SetAccess = private, Hidden = true)
        vid
        col_sum_std_fs
        row_sum_std_fs
        col_sum_lstd_sfd
        row_sum_lstd_sfd
        max_tau_fitted
    end
    
    methods
        
        %% constructor
        function obj = DDM_Analysis(filename, filepath, frameboundaries)
            % Class constructor. It loads and initializes a
            % video file, creating the initial data structure. It also
            % calculates the std of pixel intensity through time
            
            dummy_string = strsplit(filename, filesep);
            filename = dummy_string{end};
            
            if nargin < 3
                frameboundaries = [];
            end
            
            if nargin < 2 || isempty(filepath)
                if numel(dummy_string) == 1
                    filepath = pwd;
                else
                    filepath = fullfile(dummy_string{1:end-1});
                end
            end
            
            if license('test', 'Distrib_Computing_Toolbox')
                if isempty(gcp('nocreate'))
                    hcp = parpool; %starts parpool
                else
                    hcp = gcp; % uses existing parpool
                end
                hcp.IdleTimeout = 120;
            end % if license
            %% load file
            
            obj.Filename = filename;
            obj.Filepath = filepath;
            
            if isempty(frameboundaries)
                firstframe = [];
                lastframe = [];
            end
            
            fprintf('\nLoading file... ');
            
            % initialise video reader only
            obj.vid = select_video_reader(...
                fullfile(obj.Filepath, obj.Filename),...
                firstframe, lastframe);
            
            % copy properties from video reader to DDM class
            props = {'FirstFrame', 'LastFrame', 'NumberOfFrames',...
                'Height', 'Width', 'TotalNumberOfFrames',...
                'FrameRate', 'Magnification', 'px2mum'};
            for prop_ = props
                prop = prop_{1};
                obj.(prop) = obj.vid.(prop);
            end % for

            
            % read video
            % I found the global variable to be the fastest way to allow
            % all functions to access fs without it being a property
            % maybe I can write a custom 'save' command that will delete
            % the property 'fs' before saving the class on disk
            global fs;
            obj.load_movie; % writes it into fs
%             fs = obj.vid.read;
            
            % set magnification (if not read from metadata e.g. in OME)
            if isempty(obj.px2mum)
                obj.set_magnification;
            end % if

%             if isempty(strfind(filename,'.movie'))
%                 disp(filename)
%                 [fs, vi] = avi2greyscaleframestack(fullfile(obj.Filepath,obj.Filename));
%                 obj.Height      = vi.Height;
%                 obj.Width       = vi.Width;
%                 obj.FrameRate   = vi.FrameRate;
%                 if isempty(frameboundaries)
%                     obj.LastFrame   = vi.NumberOfFrames;
%                 else
%                     obj.LastFrame = min(vi.NumberOfFrames,max(frameboundaries));
%                 end
%                 obj.TotalNumberOfFrames = vi.NumberOfFrames;
%                 obj.NumberOfFrames = obj.LastFrame-obj.FirstFrame+1;
%                 obj.IsMovieCorrupted = false;
%                 fs = fs(:,:,obj.FirstFrame:obj.LastFrame);
%                 
%             elseif strfind(filename,'.movie')
%                 
%                 if ~which('moviereader')
%                     error('Cannot open .movie videos');
%                 else
%                     mo = moviereader(fullfile(obj.Filepath,obj.Filename));
%                     obj.Height      = double(mo.Height);
%                     obj.Width       = double(mo.Width);
%                     if isempty(frameboundaries)
%                         obj.LastFrame   = mo.NumberOfFrames;
%                     else
%                         obj.LastFrame = min(mo.NumberOfFrames,max(frameboundaries));
%                     end
%                     obj.TotalNumberOfFrames = mo.NumberOfFrames;
%                     
%                     obj.load_movie;
%                     %                     obj.FrameRate   = mo.FrameRate;
%                     
%                     if obj.IsMovieCorrupted    %warn user that the movie was corrupted
%                         cprintf('*[1 .3 0]',['Movie file was corrupted, only ',num2str(obj.NumberOfFrames),' frames out of ',num2str(obj.TotalNumberOfFrames),' were read. ']);
%                     end
%                 end
%                 
%             end
            cprintf('*[0 .5 0]','Done!')
            
            %% calculating standard deviation of video  % deprecated as a motion detection algorithm but let's keep it
            
            fprintf('\nCalculating standard deviation of pixel intensity through time... ');
%             std_fs = zeros(obj.Height,obj.Width,'single');
            std_fs = zeros(obj.Height,obj.Width,'double');
            parfor i=1:obj.Height    %doing the std of the entire video is too RAM expensive, so for loop on the rows
%                 std_fs(i,:) = squeeze( std( single( fs(i,:,:) ) ,1,3) );
                std_fs(i,:) = squeeze( std( double( fs(i,:,:) ) ,0,3) );
            end
            
            obj.col_sum_std_fs = sum(std_fs,1); %sum of each column of std (vertically), it's a row vector
            obj.row_sum_std_fs = sum(std_fs,2); %sum of each row of std (horizontally), it's a column vector
            
            obj.std_fs  = std_fs;   %had to do this 'cause parfor doesn't like to operate on obj.fs
            
            cprintf('*[0 .5 0]','Done!')
            
            
            %% call new motion detection algorithm
            
            obj.motion_detection;
            
            obj.col_sum_lstd_sfd = sum(obj.lstd_sfd,1);     % sum of each column of lstd_sfd (vertically), it's a row vector
            obj.row_sum_lstd_sfd = sum(obj.lstd_sfd,2);     % sum of each column of lstd_sfd (horizontally), it's a column vector
            
            
            %% initialising some global parameters
            obj.max_tau = floor(obj.NumberOfFrames/2);
            obj.max_tau_fitted = floor(obj.max_tau/2); %default, maybe change in future
            obj.set_temperature;
        end
        
        %% method that launches DDM on boxsizes specified by input vector
        function obj = VariableBoxSize_Analysis(obj, boxsizevect)
            
            if ~exist('fs','var')
                fprintf('\nAccessing global variable fs... ');
                global fs;
                if isempty(fs)
                    fprintf('fs was empty, loading the movie now... ');
                    obj.load_movie;
%                     obj.vid.read;
                end
                cprintf('*[0 .5 0]','Done!')
            end
            
            fprintf('\nStarting BoxSize Analysis... ');
            %% vector input check
            fprintf('\n\tChecking input... ');
            if isvector(boxsizevect)
                boxsizevect = boxsizevect(:); %if vector put it in column
            end
            if any(boxsizevect > obj.Width | boxsizevect > obj.Height)
                boxsizevect(boxsizevect > obj.Width | boxsizevect > obj.Height) = min(obj.Width, obj.Height);
                cprintf('*[1 .3 0]','At least one BoxSize is bigger than of the dimensions of the field of view. Shrinking it now.');
            end
            if any(mod(boxsizevect,2)) %at least one of the boxsize is odd
                boxsizevect = boxsizevect - mod(boxsizevect,2); %reduce of 1 the odd entries
                cprintf('*[1 .3 0]','At least one BoxSize is odd and this can''t be, reducing it of a unit.');
            end
            
            boxsizevect(boxsizevect == 0)=[]; %remove zero entries
            
            boxsizevect = unique(boxsizevect);  %throw away doubles, also sorts from small to big
            boxsizevect = flipud(boxsizevect); %boxsizevect is a column now, so this inverts the order
            if ~isempty(obj.Results)
                if ~isempty(intersect(boxsizevect,[obj.Results(:).BoxSize]))
                    alreadythere = intersect(boxsizevect,[obj.Results(:).BoxSize]);
                    alreadythere = alreadythere(:)'; %force to be row
                    cprintf('*[1 .3 0]',['\n\tResults for BoxSizes = ',num2str(alreadythere),' existed already, not going to overwrite it. ']);
                    fprintf('\n\t');
                    if numel(alreadythere) == numel(boxsizevect)
                        cprintf('*[1 .3 0]','\nThere''s already data for all these box sizes, skipping this video. ');
                        return
                    end
                end
                boxsizevect = setdiff(boxsizevect,[obj.Results(:).BoxSize]); %check that boxsize wasn't there already (not to overwrite results)
            end
            N_boxsizes = numel(boxsizevect);
            
            results_entry = struct('BoxSize',[],...
                'row_offset',[],...
                'col_offset',[],...
                'ind_good_boxes',[],...
                'Box',[],...
                'MedianFrequencyVec',[],...
                'AverageDamping_vs_q',[],...
                'AverageAmplitude_vs_q',[],...
                'qVec',[],...
                'MedianDamping_vs_q',[],...
                'MedianAmplitude_vs_q',[]);     %what's inside each entry of the results structure array
            
            N_existing_entries = numel(obj.Results);
            obj.Results = vertcat( obj.Results,repmat(results_entry,size(boxsizevect)));         %initialise the structure array
            clear results_entry
            
            cprintf('*[0 .5 0]','Done!')
            fprintf(['\n\tDDM will be run on Boxes of BoxSize ',num2str(boxsizevect(:)')]); % (:)' sintax to get for sure row vector
            
            %% assign BoxSize to results and call function that sets up DDM on boxes
            
            for i = 1:N_boxsizes
                j = N_existing_entries+i;
                BoxSize = boxsizevect(i);
                fprintf(['\n\tSetting up data structure for Analysis with BoxSize = ',num2str(BoxSize)]);
                obj.Results(j).BoxSize = BoxSize;
                obj.setup_DDM_on_Boxes(BoxSize);
                retrieved_max_mode_fitted = numel(obj.Results(j).Box(1).Frequency);
                obj.Results(j).qVec = 2*pi*(1:retrieved_max_mode_fitted)./BoxSize; %good for plotting, have to retrieve up until which mode I fitted the Iqtau because it is not a property of the class
                fprintf('\n\tDone');
            end
            
            % now actually run multiDDM on the boxsizes we asked
            if ~isempty(boxsizevect)
                fprintf('\n\tRunning multiDDM algorithm on all boxes... ');
                ind_to_run = N_existing_entries+1 : N_existing_entries+numel(boxsizevect);
                obj.Results(ind_to_run) = ...
                    multiDDM_core(fs, obj.N_couple_frames_to_average, obj.Results(ind_to_run));
                cprintf('*[0 .5 0]','Done!');
            end
            
            % and now fit the analysed boxes
            for i = 1:N_boxsizes
                BoxSize = boxsizevect(i);
                fprintf(['\n\tStarting fit of boxes of BoxSize = ',num2str(BoxSize)]);
                obj.fit_Boxes(BoxSize);
                fprintf('\n\tDone');
            end
            
        end
        
        
        %% method that runs DDM on all possible boxes of specified boxsize and updates correspondent results entry
        function obj = setup_DDM_on_Boxes(obj,BoxSize)
%             global fs;
            
            warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');%turn off annoying warning
            
            %% input check (function may be called directly by user, but it shouldn't)
            
            if ~isscalar(BoxSize), error('BoxSize must be a scalar - to use an array, call VariableBoxSize_Analysis instead!!');end
            if BoxSize > obj.Width || BoxSize > obj.Height
                BoxSize = min(obj.Width, obj.Height);
                cprintf('*[1 .3 0]',['Impossible to run the analysis with the desired BoxSize, as it is bigger than at least one of the dimensions of the field of view. The analysis will be run on the maximum possible BoxSize, i.e. BoxSize = ',num2str(BoxSize)]);
            end
            if mod(BoxSize,2) %if odd BoxSize reduce of 1
                cprintf('*[1 .3 0]','BoxSize can''t be odd, reducing it of a unit. ');
                BoxSize = BoxSize-1;
            end
            if BoxSize == 0
                error('BoxSize can''t be 0.');
            end
            
            %% now find where to put the results
            
            if isempty(obj.Results)
                iir = 1;
                fprintf(['\n\t\tCreating a result entry for the desired BoxSize... Results will be stored in entry ',num2str(iir)]);
                obj.Results(iir).BoxSize = BoxSize;
            else
                iir = find([obj.Results(:).BoxSize] == BoxSize); %look for an entry of results with desired BoxSize
                if numel(iir) > 1
                    error('Something went wrong in input checking. There are two results entries with the same BoxSize')
                end
                if isempty(iir)
                    iir = numel(obj.Results)+1;
                    fprintf(['\n\t\tCreating a result entry for the desired BoxSize... Results will be stored in entry ',num2str(iir)]);
                    obj.Results(iir).BoxSize = BoxSize;
                    
                end
            end
            
            %% prediction on sizes
            
            N_Boxes_row = floor(obj.Height / BoxSize);   %number of boxes that fit (vertically) in the field of view given BoxSize
            N_Boxes_col = floor(obj.Width  / BoxSize);    %number of boxes that fit (horizontally) in the field of view given BoxSize
            N_Boxes = N_Boxes_row * N_Boxes_col;        %number of boxes that fit (in total) in the field of view given BoxSize
            
            max_mode = BoxSize/2;           %number of rows of Iqtau
            
            max_mode_fitted = floor(max_mode/2);    %maximum mode (higher q) fitted
            
            
            %% initialisation of result structure
            
            Box_struct = struct(...
                'std_fs',    zeros(BoxSize,'single'),...
                'Iqtau',     zeros(max_mode, obj.max_tau, 'single'),...
                'Frequency', zeros(max_mode_fitted,1),...
                'Amplitude', zeros(max_mode_fitted,1),...
                'Damping',   zeros(max_mode_fitted,1),...
                'Offset',    zeros(max_mode_fitted,1),...
                'GOF',       zeros(max_mode_fitted,1));
            obj.Results(iir).Box = repmat(Box_struct, [N_Boxes, 1]);
            clear Box_struct;
            
            %% choosing offsets to maximise output
            % all the boxes touch each other, but if the sizes of the video are not
            % multiples of BoxSize there is a bit of wiggle room to play with to
            % maximise the signal we're analysing. this section just uses the standard
            % deviation of the video to detect the most active region
            
            fprintf('\n\t\tChoosing ROI placement... ');
            
            col_span = N_Boxes_col * BoxSize;     %Width of total DDM-analysed region (multiple of BoxSize)
            row_span = N_Boxes_row * BoxSize;    %Height of total DDM-analysed region (multiple of BoxSize)
            
            %find where to place (in the horizontal direction) the [row_span, col_span] region to be DDM-analysed
            for i = obj.Width - col_span + 1 : -1 : 1  %sneaky allocation
%                 dummy(i) = sum(obj.col_sum_std_fs(i:i+col_span-1));
                dummy(i) = sum(obj.col_sum_lstd_sfd(i:i+col_span-1));
            end
            [~,col_offset] = max(dummy);
            clear dummy
            %find where to place (in the vertical direction) the [row_span, col_span] region to be DDM-analysed
            for i = obj.Height - row_span + 1 : -1 : 1  %sneaky allocation
%                 dummy(i) = sum(obj.row_sum_std_fs(i:i+row_span-1));
                dummy(i) = sum(obj.row_sum_lstd_sfd(i:i+row_span-1));
            end
            [~,row_offset] = max(dummy);
            clear dummy
            cprintf('*[0 .5 0]','Done!')
            
            %% select good boxes based on thresholding of standard deviation
            
%             ind_good_boxes =
%             find_good_boxes(obj.std_fs(row_offset:row_offset+row_span-1,
%             col_offset:col_offset+col_span-1), BoxSize); %old function
            ind_good_boxes = obj.find_boxes_with_motion(row_offset, col_offset, BoxSize);
            
            %% saving in each Box the relative bit of std_fs
            
%             fprintf('\n\t\tSaving portion of standard deviation of pixel intensity in time relative to each Box... ');
            for i = 1:N_Boxes_row
                for j = 1:N_Boxes_col
                    
                    ii = sub2ind([N_Boxes_row, N_Boxes_col], i, j);
                    obj.Results(iir).Box(ii).std_fs = [];
                    obj.Results(iir).Box(ii).lstd_sfd = [];
%                     obj.Results(iir).Box(ii).std_fs = obj.std_fs(row_offset + (i-1)*BoxSize : row_offset + (i)*BoxSize -1 ,...
%                         col_offset + (j-1)*BoxSize : col_offset + (j)*BoxSize -1);
%                     obj.Results(iir).Box(ii).lstd_sfd = obj.lstd_sfd(row_offset + (i-1)*BoxSize : row_offset + (i)*BoxSize -1 ,...
%                         col_offset + (j-1)*BoxSize : col_offset + (j)*BoxSize -1);
                end
            end
%             cprintf('*[0 .5 0]','Done!')
            
            %% running DDM on each region of the video, storing Iqtaus in relative Box
            % not done here anymore, see end of Variable_BoxSize_analysis
            %{
            fprintf('\n\t\tRunning DDM algorithm on each Box... ');
            fprintf('\n\t\t\tAnalysing Box ');
            for j = 1:N_Boxes_col
                for i = 1:N_Boxes_row
                    
                    ii = sub2ind([N_Boxes_row, N_Boxes_col], i, j);
                    fprintf('%.6d / %.6d',ii,N_Boxes);
                    obj.Results(iir).Box(ii).Iqtau = DDM_core(fs(row_offset + (i-1)*BoxSize : row_offset + (i)*BoxSize -1 ,...
                        col_offset + (j-1)*BoxSize : col_offset + (j)*BoxSize -1,:),obj.N_couple_frames_to_average);
                    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
                end
            end
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % deletes "Analysing Box"
            
            cprintf('*[0 .5 0]','Done!');
            %}
            
            %% running fit_Iqtau_ourcilia on each Box.Iqtau
            % not done here anymore, see end of Variable_BoxSize_analysis
            %{
            fprintf('\n\t\tFitting Iqtau matrix for each Box... ');
            fprintf('\n\t\t\tFitting Box ');
            for ii = 1:N_Boxes
                
                fprintf('%.6d / %.6d',ii,N_Boxes);
                [obj.Results(iir).Box(ii).Frequency, obj.Results(iir).Box(ii).Damping, obj.Results(iir).Box(ii).Amplitude, obj.Results(iir).Box(ii).Offset, obj.Results(iir).Box(ii).GOF] = fit_Iqtau_ourcilia(obj.Results(iir).Box(ii).Iqtau, max_mode_fitted, obj.max_tau_fitted);
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
            end
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % deletes "fitting box"
            
            cprintf('*[0 .5 0]','Done!');
            %}
            
            %% saving parameters in results
            
            obj.Results(iir).row_offset = row_offset;
            obj.Results(iir).col_offset = col_offset;
            obj.Results(iir).ind_good_boxes = ind_good_boxes;
            
        end
        
        %% function that runs fits on the boxes
        function obj = fit_Boxes(obj, BoxSize)
            
            % find where to put the results
            iir = find([obj.Results(:).BoxSize] == BoxSize); %look for an entry of results with desired BoxSize
            if numel(iir) > 1
                error('Something went wrong in input checking. There are two results entries with the same BoxSize')
            end
            if isempty(iir)
                error('The boxes of the BoxSize you wish to fit haven''t been analysed by multiDDM yet');
            end

            % running fit_Iqtau_ourcilia on each Box.Iqtau
            N_Boxes = numel(obj.Results(iir).Box);
            max_mode_fitted = floor(BoxSize/4);
            
            fprintf('\n\t\tFitting Iqtau matrix for each Box... ');
            fprintf('\n\t\t\tFitting Box ');
            for ii = 1:N_Boxes
                
                fprintf('%.6d / %.6d',ii,N_Boxes);
                [freq, damp, ampl, offs, gof] = fit_Iqtau_ourcilia(obj.Results(iir).Box(ii).Iqtau,...
                    max_mode_fitted, obj.max_tau_fitted);
                obj.Results(iir).Box(ii).Frequency = freq;
                obj.Results(iir).Box(ii).Damping = damp;
                obj.Results(iir).Box(ii).Amplitude = ampl;
                obj.Results(iir).Box(ii).Offset = offs;
                obj.Results(iir).Box(ii).GOF = gof;
                
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
            end
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % deletes "fitting box"
            
            cprintf('*[0 .5 0]','Done!');
            
        end %function
        
        %% interactive creation of kymograph
        function obj = kymograph(obj,N_kymographs)
            
            if nargin==1
                N_kymographs = 1;
            end
            
            global fs;
            if size(fs) ~= horzcat(obj.Height, obj.Width, obj.NumberOfFrames)
                cprintf('*[1 .3 0]','It looks like the video in the memory is not the same as specified by the instance of the class. Loading the proper video...');
%                 fs = obj.vid.read;
                obj.load_movie;
            end
            
            fprintf('\nAll previous kimographs for this video will be deleted!\n');
            
            obj.Kymograph = repmat(struct('ind',[], 'kymo',[]),[N_kymographs,1]);
            
            hf = figure(11001);
            imshow(fs(:,:,1),[]);
            
            for ki = 1:N_kymographs
                [col_kind, row_kind, ffprofile] = improfile;
                col_kind = round(col_kind);
                row_kind = round(row_kind);
                obj.Kymograph(ki).ind = sub2ind([obj.Height, obj.Width], row_kind, col_kind);
                obj.Kymograph(ki).kymo = zeros(numel(col_kind),obj.NumberOfFrames,'double');
                
            end
            
            for ki = 1:N_kymographs
                for t = 1:obj.NumberOfFrames
                    im = fs(:,:,t);
                    obj.Kymograph(ki).kymo(:,t) = im(obj.Kymograph(ki).ind);
                end
            end
            
            
            close(hf);
            
            
        end
        
        %% method that loads the movie (useful if fs got cleared)
        function obj = load_movie(obj)
            %load_movie wrapper for internal reader's read method
            global fs;
            fs = obj.vid.read;
        end % function
        
%         function [obj, movie_correctly_read] = load_movie(obj)
%             
%             movie_correctly_read = false;   %control variable for failsafe movie reading
%             
%             global fs;
%             if isprop(obj,'Filepath')
%                 fullpath = fullfile(obj.Filepath,obj.Filename);
%             else
%                 fullpath = obj.Filename;
%             end
%             
%             [~,~,ext] = fileparts(fullpath);
%             
%             if strcmp(ext,'.movie')
%                 try
%                     mo = moviereader(fullpath);
%                 catch
%                     cprintf('Red','The movie is not present at the path specified ');
%                     return
%                 end
%                 
%                 obj.IsMovieCorrupted = false;
%                 
%                 while movie_correctly_read == false %try to read the movie until it works
%                     try
%                         fs = mo.read([obj.FirstFrame obj.LastFrame]);  %try to read the movie
%                         movie_correctly_read = true;        %if it works put flag to 1
%                         
%                     catch err                               %if it doesn't work, then
%                         if strcmp(err.identifier,'moviereader:WrongMagicAtFrame')
%                             obj.LastFrame = str2double(err.message)-1;  %try to read again, ditching the last frame (I guess I could use a "finding" algorithm here)
%                         else
%                             obj.LastFrame = obj.LastFrame - 1;
%                         end
%                         obj.IsMovieCorrupted = true;        %set the flag that the movie was corrupted
%                         
%                     end                                     %and try again
%                     
%                 end
%                 
%                 
%                 fs(:,:,end) = []; %ditches the last frame anyway, which may still be corrupted
%                 obj.LastFrame = obj.LastFrame - 1;
%                 
%                 if obj.IsMovieCorrupted
%                     cprintf('Red','The movie is corrupted ')
%                 end
%                 
%                 
%                 
%                 obj.FrameRate   = mo.FrameRate;
%                 obj.NumberOfFrames = obj.LastFrame - obj.FirstFrame + 1;
%                 
%                 
%                 
%             else
%                 fs = avi2greyscaleframestack(fullpath);
%                 fs = fs(:,:,obj.FirstFrame:obj.LastFrame);
%             end
%             
%             if obj.NumberOfFrames ~= size(fs,3)
%                 error('mistake here')
%             end
%             
%         end
        
        %% parse filename for magnification string
        function obj = set_magnification(obj)
            
            [px2mum_, magstr] = magnification_to_micronsperpixel(...
                fullfile(obj.Filepath, obj.Filename),...
                'parameters\calibrated_magnifications.json');
            
            if px2mum_ < 0
                warning('No valid micron/pixel. Will continue anyway')
            else
                obj.px2mum = px2mum_;
                obj.Magnification = magstr;
            end %if
            
        end
        
        %% function that gathers the results (does means etc)
        function obj = gather_results(obj, flag_recalculate_good_boxes)
            
            % input checking
            
            if nargin < 2 || isempty(flag_recalculate_good_boxes)
                flag_recalculate_good_boxes = false;
            end
            
            
            for ii = 1:numel(obj.Results) %create a vector with median frequency of each box: as long as the number of boxes for each boxsize
                
                if ~isfield(obj.Results(ii),'ind_good_boxes') ||...
                        (isfield(obj.Results(ii),'ind_good_boxes') && isempty(obj.Results(ii).ind_good_boxes) ) || ...
                        flag_recalculate_good_boxes
%                     row_span = floor(obj.Height / obj.Results(ii).BoxSize) * obj.Results(ii).BoxSize;   %number of boxes that fit (vertically) in the field of view given BoxSize
%                     col_span = floor(obj.Width  / obj.Results(ii).BoxSize) * obj.Results(ii).BoxSize;
%                     ind_good_boxes = find_good_boxes(obj.std_fs(obj.Results(ii).row_offset:obj.Results(ii).row_offset+row_span-1,obj.Results(ii).col_offset:obj.Results(ii).col_offset+col_span-1), ...
%                         obj.Results(ii).BoxSize, flag_ditch_imadjust);
                    ind_good_boxes = obj.find_boxes_with_motion(obj.Results(ii).row_offset,...
                        obj.Results(ii).col_offset, ...
                        obj.Results(ii).BoxSize);
                    obj.Results(ii).ind_good_boxes = ind_good_boxes;
                end
                
                temp_Frequency_mat = horzcat(obj.Results(ii).Box(obj.Results(ii).ind_good_boxes).Frequency).*obj.FrameRate; % exclude bad boxes (threshold of std_fs) first
                temp_Damping_mat = horzcat(obj.Results(ii).Box(obj.Results(ii).ind_good_boxes).Damping).*obj.FrameRate;
                temp_Amplitude_mat = horzcat(obj.Results(ii).Box(obj.Results(ii).ind_good_boxes).Amplitude);
                
                bad_fits = false(size(temp_Frequency_mat)); %controllo per buttar via roba
                
                bad_fits(temp_Frequency_mat <=   1)  = true;    %exclude boxes where frequency too low (aka sth went wrong in the fit or the sample is still)
                bad_fits(temp_Frequency_mat >=  30)  = true;    %exclude boxes where frequency too high (aka sth went wrong in the fit)
                bad_fits(temp_Damping_mat   ==   0)  = true;    %exclude boxes where there is no damping (aka sth went wrong in the fit)
                bad_fits(temp_Amplitude_mat <= eps)  = true;    %exclude boxes where signal too low
                
                temp_Frequency_mat(bad_fits) = NaN;
                temp_Amplitude_mat(bad_fits) = NaN;
                temp_Damping_mat(bad_fits) = NaN;
                
                %                 obj.Results(ii).MedianFrequencyVec = arrayfun(@(bc)nanmedian(obj.Results(ii).Box(bc).Frequency), 1:numel(obj.Results(ii).Box)); % bc = box_counter
                %                 obj.Results(ii).AverageDamping_vs_q = nanmean( horzcat(obj.Results(ii).Box(:).Damping) ,2);
                %                 obj.Results(ii).AverageAmplitude_vs_q = nanmean( horzcat(obj.Results(ii).Box(:).Amplitude), 2);
                
                obj.Results(ii).MedianFrequencyVec = prctile(temp_Frequency_mat, 50); % bc = box_counter %% IN HERTZ NOW!!!!
                obj.Results(ii).AverageDamping_vs_q = trimmean(temp_Damping_mat, 20, 2);   %% IN HERTZ NOW!!!!
                obj.Results(ii).AverageAmplitude_vs_q = trimmean(temp_Amplitude_mat, 10, 2);
                
                obj.Results(ii).MedianDamping_vs_q = prctile(temp_Damping_mat, 50, 2);   %% IN HERTZ NOW!!!!
                obj.Results(ii).MedianAmplitude_vs_q = prctile(temp_Amplitude_mat, 50, 2);
                
                if ~isfield(obj.Results,'qVec') %should enter here only for ii=1
                    retrieved_max_mode_fitted = numel(obj.Results(ii).AverageDamping_vs_q);
                    obj.Results(ii).qVec = 2*pi*(1:retrieved_max_mode_fitted)./obj.Results(ii).BoxSize;
                else
                    if isempty(obj.Results(ii).qVec) %when you create qVec for ii=1 it is created empty of all other entries of Results, so it will exist but be empty
                        retrieved_max_mode_fitted = numel(obj.Results(ii).AverageDamping_vs_q);
                        obj.Results(ii).qVec = 2*pi*(1:retrieved_max_mode_fitted)./obj.Results(ii).BoxSize;
                    end
                end
            end
            
            
            if ~isfield(obj.SAVAlike,'ind_good_bins') && ~isempty(obj.SAVAlike)
                bsz = mean(floor([obj.Height, obj.Width]./size(obj.SAVAlike.frequency_map)));
%                 obj.SAVAlike.ind_good_bins = find_good_boxes(obj.std_fs,bsz);
                obj.SAVAlike.ind_good_bins = obj.find_boxes_with_motion(1,1,bsz);
                ind = obj.SAVAlike.frequency_map(:) > 1 & obj.SAVAlike.frequency_map(:) < 30; %physical constraints
                ind = ind & obj.SAVAlike.ind_good_bins(:);          % discard empty bins by thresholding the std
                obj.SAVAlike.mean_frequency = mean(obj.SAVAlike.frequency_map(ind));
                obj.SAVAlike.std_frequency = std(obj.SAVAlike.frequency_map(ind));
                obj.SAVAlike.err_frequency = obj.SAVAlike.std_frequency/sqrt(sum(ind)); %error of the mean
                obj.SAVAlike.median_frequency = median(obj.SAVAlike.frequency_map(obj.SAVAlike.ind_good_bins));   %shouldn't be too different than the constrained mean a few lines above
            end
            
        end
        
        %% frequency_map with SAVAlike software
        function obj = SAVAlike_CBF_measurement(obj, flag_ROI)
            if nargin < 2, flag_ROI = false; end
            bsz = 4; %size of binning box
            global fs;
            
            fprintf('\nStarting SAVA-like CBF measurement...');
            
            
            rect = [1 1 obj.Width obj.Height];
            
            if flag_ROI %select ROI
                figure(11002);
                imagesc(obj.std_fs); colormap jet;
                
                rect = round(getrect(gcf));
            end
            
            % resize fs so that frame sizes are a multiple of bsz 
            rect(3) = bsz*ceil(rect(3)/bsz);
            rect(4) = bsz*ceil(rect(4)/bsz);

            fs = fs(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1,:);
            
            
            ROI_height = rect(4);
            ROI_width  = rect(3);
            
            % binning
            fprintf('\n\tBinning...');
            binningmap = col2im(repmat(1:floor(ROI_width*ROI_height/bsz^2),bsz*bsz,1), bsz.*[1 1], bsz*[floor(ROI_height/bsz) floor(ROI_width/bsz)],'distinct');
            binnedfs = zeros(floor( [ROI_height/bsz ROI_width/bsz obj.NumberOfFrames]), 'single');
            binna = @(i) reshape( 1./bsz^2.*accumarray(binningmap(:),single(reshape(fs(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1,i),ROI_width*ROI_height,1))), [floor(ROI_height/bsz), floor(ROI_width/bsz)]);
            for i=1:obj.NumberOfFrames, binnedfs(:,:,i) = binna(i); end
            cprintf('*[0 .5 0]','Done!')
            
            % initialisation of variable for spectrum calculations
            bheight = floor(size(binningmap,1)/bsz);
            bwidth  = floor(size(binningmap,2)/bsz);
            frequency_map = zeros(bheight, bwidth,'double');
            
            % calculation of CBF
            fprintf('\n\tMeasuring CBF...');
            binnedfs = reshape(binnedfs,[bheight*bwidth obj.NumberOfFrames]);
            fy = abs(fft(binnedfs,[],2));
            spectrum1 = fy(:,1:ceil(obj.NumberOfFrames/2)); %next three line from algorithm to detect peaks in signal with harmonics
            spectrum2 = fy(:,1:2:obj.NumberOfFrames);
            
            parfor ii=1:bheight*bwidth
                spectrum1(ii,:)=smooth(spectrum1(ii,:));
                spectrum2(ii,:)=smooth(spectrum2(ii,:));
            end
            
            product_spectrum = spectrum1.*spectrum2;
            
            %             pks_map = zeros(size(frequency_map),'double');
            fr = obj.FrameRate; %for parfor reasons
            nf = obj.NumberOfFrames;
            parfor ii=1:bheight*bwidth
                
                [pks,loc] = findpeaks( double(product_spectrum(ii,:)) );
                %     [~,loc,~,proms] = findpeaks( double(spectrum1(ii,:)) ); %good only for very clean spectra
                %                 [pkmax,ind] = max(pks);
                [~,ind] = max(pks);
                %     [loc, pks] = peakseek(product_spectrum(ii,:));
                %     [~,ind] = max(pks);
                if ~isempty(ind) %to prevent error when there-s a saturated pixel
                    frequency_map(ii) = (loc(ind)+1)/(nf).* fr; %converting to Hz (? to be checked)
                end
                %                 pks_map(ii) = pkmax;
            end
            cprintf('*[0 .5 0]','Done!')
            
            % results storage
            obj.SAVAlike.frequency_map = frequency_map;
            obj.SAVAlike.ind_good_bins = find_good_boxes(obj.std_fs,bsz);
            ind = frequency_map(:) > 1 & frequency_map(:) < 30; %physical constraints
            ind = ind & obj.SAVAlike.ind_good_bins(:);          % discard empty bins by thresholding the std
            obj.SAVAlike.mean_frequency = mean(frequency_map(ind));
            obj.SAVAlike.std_frequency = std(frequency_map(ind));
            obj.SAVAlike.err_frequency = obj.SAVAlike.std_frequency/sqrt(sum(ind)); %error of the mean
            obj.SAVAlike.median_frequency = median(frequency_map(:));   %shouldn't be too different than the constrained mean a few lines above
            cprintf('*[0 .5 0]','\nDone!')
            
        end
        
        %% read temperature file
        function obj = set_temperature(obj, channel)
            %set_temperature read temperature from companion case (only if
            %video was in .movie format)
            [~,~,movie_ext] = fileparts(obj.Filename);
            if ~strcmp(movie_ext, '.movie')
                return
            end % if
            
            if nargin < 2
                channel = 0; % select TC1
            end
            channel = channel + 2; %to take into account the time vector column and C notation
            
            TemperatureFilename = strrep(obj.Filename, movie_ext, '.temperature');
            if isprop(obj,'Filepath')
                TemperatureFilename = fullfile(obj.Filepath,TemperatureFilename);
            end
            if isempty(dir(TemperatureFilename))
                fprintf('\n');
                cprintf('*[1 .3 0]',regexprep(['Temperature file not found at ',TemperatureFilename],'\','\\\\'));
            else
                try
                    dummy = importdata(TemperatureFilename);
                    if ~isempty(dummy)
                        %assuming the temperature time vector is ~ as long as true
                        %time vector
                        dt = obj.NumberOfFrames/obj.FrameRate/size(dummy.data,1);
                        obj.Temperature.TimeVec = dt * (dummy.data(:,1) - dummy.data(1,1)); %time vector (in the temperature file is in units of 0.1s and with god-knows-what offset)
                        obj.Temperature.TempVec = dummy.data(:,channel);
                        obj.Temperature.MeanTemp = mean(obj.Temperature.TempVec);
                    else
                        cprintf('*[1 .3 0]','Temperature file empty.');
                    end
                catch
                    cprintf('*[1 .3 0]','Reading Temperature file failed.');
                end
            end
        end
        
        
        
        %% better movement finder       
        function obj = motion_detection(obj, flag_legacy)
            % starting from the first second of video, find areas with
            % cilia beating
            
            if nargin < 2 || isempty(flag_legacy)
                flag_legacy = false;
            end
            
            sfwsz = 11;                 % size of the stdfiltering window
                
            nbins_greylevels = 512;     % number of bins for grey level histogram

            if isempty(obj.lstd_sfd)
                
                
                % check if video is loaded, else load the first 1 second
                if ~exist('fs','var')
                    fprintf('\nAccessing global variable fs... ');
                    global fs;
                    if isempty(fs)
                        fprintf('fs was empty, reading the first 1s of the movie now... ');
                    
                        % read the first second
                        try
                            fs = obj.vid.read([1 round(obj.FrameRate)]);
                        catch EE
                            disp('Error trying to read:')
                            disp(fullfile(obj.Filepath, obj.Filename))
                            error(EE.message)
                        end
                    end % if isempty
                    
                end %if ~exist
                
                % initialise stack
                nf = min(size(fs,3), round(obj.FrameRate)) -2; % 1 s worth, or all video
                sfd = zeros([obj.Height, obj.Width, nf], 'single'); %stdfiltered diffs
                for i = size(sfd,3):-1:1
                    sfd(:,:,i) = single(stdfilt( medfilt2( diff(single(fs(:,:,[i,i+2])),1,3) ,'symmetric') , ones(sfwsz)));
                end %for
                
                % now take std in time
                std_sfd = std(sfd,0,3);
                std_sfd(std_sfd == 0) = min(std_sfd(std_sfd>0)); % prevents an error if the video was saturated
                obj.lstd_sfd = mat2gray(log10( std_sfd ));       % log on graylevels, deals better with peak colours
                
                % fix the very high values at the very edges
                mm = min(obj.lstd_sfd(:));
                obj.lstd_sfd([1:ceil(sfwsz/2), end-floor(sfwsz/2):end], :) = mm;
                obj.lstd_sfd(:, [1:ceil(sfwsz/2), end-floor(sfwsz/2):end]) = mm;
                
                
            end %if no log motion map
            
            if flag_legacy
                
                cprintf('*[1 .3 0]','Using legacy motion detection')
                level = graythresh(obj.lstd_sfd(ceil(sfwsz/2)+1:end-ceil(sfwsz/2),ceil(sfwsz/2)+1:end-ceil(sfwsz/2)));
                level = 0.7 * level;
                obj.detected_motion = imopen(imbinarize(obj.lstd_sfd, level),strel('disk',9));
                
            else
                
                % fit the histogram of grey levels of lstd_sfd
                % with a gaussian, then put to minimum all the points
                % within a sigms of the main peak which is usually the
                % background
                [yhist, xhist] = histcounts(obj.lstd_sfd(:), linspace(0, 1, nbins_greylevels + 1),...
                    'normalization','probability');
                yhist = yhist(:);
                xhist = xhist(:);
                xhist = (xhist(1:end-1)+xhist(2:end))/2;
                
                % initial values for fitting
                [a0, b0] = findpeaks(smooth(yhist),xhist,...
                    'NPeaks',1,'MinPeakHeight',0.5e-3);
                c0 = 0.1;
                
                % fitting
                ft = fittype('gauss1');
                fo = fitoptions('method','nonlinearleastsquares',...
                    'startpoint',[a0, b0, c0],...
                    'Upper',[1 1 1],'lower',[0 0 0]);
                idx_to_fit = xhist<=b0+c0/sqrt(2);
                idx_to_fit(1) = false;  %because of edges, this is overly populated
                fit_out = fit(xhist(idx_to_fit),yhist(idx_to_fit),ft,fo);
                sgm = fit_out.c1/sqrt(2);
                
                % now find background, put it to min
                bg = obj.lstd_sfd <= fit_out.b1+2*sgm;
                bg = imclose(bg,strel('disk',9));
                
                
                % if this failed call legacy method
                if any(~bg(:)) % all's well
                    % binarise
                    obj.detected_motion = ~bg;
                else
                    obj.motion_detection(true);
                end
                
            end % if flag_legacy
            
        end %function
        
        
        %% find where the motion mask is true within a box
        function mask = find_boxes_with_motion(obj, row_offset, col_offset, bsz)
            %find_boxes_with_motion checks which ones of the bsz-by-bsz
            %boxes is over a region detected as true by detect_motion
            
            % run motion detection if needed
            if isempty(obj.detected_motion) || isempty(obj.lstd_sfd)
                obj.motion_detection;
            end
            
            row_span = floor(obj.Height / bsz) * bsz;   %number of boxes that fit (vertically) in the field of view given BoxSize
            col_span = floor(obj.Width  / bsz) * bsz;

            % chop the edges of detected motion off according to row_offset
            % and col_offset
            BW = obj.detected_motion(row_offset:row_offset+row_span-1,col_offset:col_offset+col_span-1);
            
            % prepare binning map to go from bigmask to mask, this is
            % still common to 2 3 4 methods
            w = size(BW,2);
            h = size(BW,1);
            binningmap = col2im(repmat(1:floor(w*h/bsz^2),bsz*bsz,1),...
                bsz.*[1 1], bsz*[floor(h/bsz) floor(w/bsz)],'distinct');
            
            % take the box as long as there is 1 true pixel in the mask under it
            mask = reshape(accumarray(binningmap(:), BW(:)), h/bsz, w/bsz) > 0;
            
        end %function
        
    end %methods

    
end

%% select good boxes based on thresholding of standard deviation

function mask = find_good_boxes(IM, bsz, flag_ditch_imadjust)

if nargin < 3 || isempty(flag_ditch_imadjust)
    flag_ditch_imadjust = false;
end %if

IM = mat2gray(IM);      %converts to [0 1]
[height, width] = size(IM);
binningmap = col2im(repmat(1:floor(width*height/bsz^2),bsz*bsz,1), bsz.*[1 1], bsz*[floor(height/bsz) floor(width/bsz)],'distinct');


if flag_ditch_imadjust

    % highlight local roughness
    fIM = stdfilt(IM, ones(5));
    
    % smooth it
    sIM = mat2gray(imgaussfilt(fIM, 5));
    saIM = imadjust(sIM, stretchlim(sIM, [0.01, 1-1e-3]));
    
    % global thresholding
    [level, ~] = graythresh(saIM);
    bigmask = imbinarize(saIM, 0.65*level);
    bigmask = bwmorph(bigmask,'clean');

    % take the box as long as there is 1 true pixel in the mask under it
    mask = reshape(accumarray(binningmap(:), bigmask(:)), height/bsz, width/bsz) > 0;
    
else
    
    bIM = reshape( accumarray(binningmap(:),IM(:)), height/bsz, width/bsz) ./ bsz.^2;
    abIM = imadjust(bIM);
    mask = im2bw(abIM,min(1,graythresh(bIM)));

end %if


if bsz < 64
    mask = bwmorph(mask,'clean'); % maybe want to think about this again
end

end




