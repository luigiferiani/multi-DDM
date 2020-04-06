%{
% Version 1.0
% © Luigi Feriani, Maya Juenet, Pietro Cicuta, 2017 (lf352@cam.ac.uk) 
% 
% DDM_Analysis.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% 
% Original work
% Feriani, L., et al., Biohpysical Journal 2017
% "Assessing the Collective Dynamics of Motile Cilia in Cultures of Human
% Airway Cells by Multiscale DDM"
% http://dx.doi.org/10.1016/j.bpj.2017.05.028
%
%}

classdef DDM_Analysis < matlab.mixin.Copyable
    %DDM_analysis class to handle and run ddm analysis on microscope videos
    %   Detailed explanation goes here
    
    properties
        Filename            %Name of the video file
        Filepath            %Name of file path
        height              %rows of video, (px)
        width               %cols of video, (px)
        TotalNumberOfFrames %N_frames as written in the movie header
        FirstFrame          %first frame to analyse
        LastFrame           %last frame to analyse
        NumberOfFrames      %number of frames to analyse (LastFrame-FirstFrame+1)
        FrameRate           %
        Magnification       %for px2mum conversions
        std_fs              %standard deviation of pixel intensity through time
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
    end
    
    properties (SetAccess = private, Hidden = true)
        col_sum_std_fs
        row_sum_std_fs
        max_tau
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
            
            if nargin < 2
                if numel(dummy_string) == 1
                    filepath = pwd;
                else
                    filepath = fullfile(dummy_string{1:end-1});
                end
            end
            
            if isempty(gcp('nocreate'))
                hcp = parpool; %starts parpool
            else
                hcp = gcp; % uses existing parpool
            end
            hcp.IdleTimeout = 120;
            %% load file
            
            obj.Filename = filename;
            obj.Filepath = filepath;
            
            if isempty(frameboundaries)
                obj.FirstFrame = 1;
            else
                obj.FirstFrame = max(1,min(frameboundaries));
            end
            
            fprintf('\nLoading file... ');
            
            global fs;
            if isempty(strfind(filename,'.movie'))
                filename
                [fs, vi] = avi2greyscaleframestack(fullfile(obj.Filepath,obj.Filename));
                obj.height      = vi.Height;
                obj.width       = vi.Width;
                obj.FrameRate   = vi.FrameRate;
                if isempty(frameboundaries)
                    obj.LastFrame   = vi.NumberOfFrames;
                else
                    obj.LastFrame = min(vi.NumberOfFrames,max(frameboundaries));
                end
                obj.TotalNumberOfFrames = vi.NumberOfFrames;
                obj.NumberOfFrames = obj.LastFrame-obj.FirstFrame+1;
                obj.IsMovieCorrupted = false;
                fs = fs(:,:,obj.FirstFrame:obj.LastFrame);
                
            elseif strfind(filename,'.movie')
                
                if ~which('moviereader')
                    error('Cannot open .movie videos');
                else
                    mo = moviereader(fullfile(obj.Filepath,obj.Filename));
                    obj.height      = double(mo.height);
                    obj.width       = double(mo.width);
                    if isempty(frameboundaries)
                        obj.LastFrame   = mo.NumberOfFrames;
                    else
                        obj.LastFrame = min(mo.NumberOfFrames,max(frameboundaries));
                    end
                    obj.TotalNumberOfFrames = mo.NumberOfFrames;
                    
                    obj.load_movie;
                    %                     obj.FrameRate   = mo.FrameRate;
                    
                    if obj.IsMovieCorrupted,    %warn user that the movie was corrupted
                        cprintf('*[1 .3 0]',['Movie file was corrupted, only ',num2str(obj.NumberOfFrames),' frames out of ',num2str(obj.TotalNumberOfFrames),' were read. ']);
                    end
                end
                
            end
            cprintf('*[0 .5 0]','Done!')
            
            %% calculating standard deviation of video
            
            fprintf('\nCalculating standard deviation of pixel intensity through time...');
%             std_fs = zeros(obj.height,obj.width,'single');
            std_fs = zeros(obj.height,obj.width,'double');
            parfor i=1:obj.height    %doing the std of the entire video is too RAM expensive, so for loop on the rows
%                 std_fs(i,:) = squeeze( std( single( fs(i,:,:) ) ,1,3) );
                std_fs(i,:) = squeeze( std( double( fs(i,:,:) ) ,1,3) );
            end
            
            obj.col_sum_std_fs = sum(std_fs,1); %sum of each column of std (vertically), it's a row vector
            obj.row_sum_std_fs = sum(std_fs,2); %sum of each row of std (horizontally), it's a column vector
            
            obj.std_fs  = std_fs;   %had to do this 'cause parfor doesn't like to operate on obj.fs
            
            cprintf('*[0 .5 0]','Done!')
            
            %% initialising some global parameters
            obj.max_tau = floor(obj.NumberOfFrames/2);
            obj.max_tau_fitted = floor(obj.max_tau/2); %default, maybe change in future
            obj.set_magnification;
            obj.set_temperature;
        end
        
        %% method that launches DDM on boxsizes specified by input vector
        function obj = VariableBoxSize_Analysis(obj,boxsizevect)
            
            fprintf('\nStarting BoxSize Analysis...');
            %% vector input check
            fprintf('\n\tChecking input...');
            if isvector(boxsizevect),
                boxsizevect = boxsizevect(:); %if vector put it in column
            end
            if any(boxsizevect > obj.width | boxsizevect > obj.height)
                boxsizevect(boxsizevect > obj.width | boxsizevect > obj.height) = min(obj.width, obj.height);
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
                if ~isempty(intersect(boxsizevect,[obj.Results(:).BoxSize])),
                    alreadythere = intersect(boxsizevect,[obj.Results(:).BoxSize]);
                    alreadythere = alreadythere(:)'; %force to be row
                    cprintf('*[1 .3 0]',['\nResults for BoxSizes = ',num2str(alreadythere),' existed already, not going to overwrite it.']);
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
            
            %% assign BoxSize to results and call function that runs DDM on boxes
            for i = 1:N_boxsizes
                j=N_existing_entries+i;
                BoxSize = boxsizevect(i);
                fprintf(['\n\tStarting Analysis with BoxSize = ',num2str(BoxSize)]);
                obj.Results(j).BoxSize = BoxSize;
                obj.DDM_on_Boxes(BoxSize);
                retrieved_max_mode_fitted = numel(obj.Results(j).Box(1).Frequency);
                obj.Results(j).qVec = 2*pi*[1:retrieved_max_mode_fitted]./BoxSize; %good for plotting, have to retrieve up until which mode I fitted the Iqtau because it is not a property of the class
                fprintf(['\n\tDone']);
            end
            
        end
        
        
        %% method that runs DDM on all possible boxes of specified boxsize and updates correspondent results entry
        function obj = DDM_on_Boxes(obj,BoxSize)
            global fs;
            
            warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');%turn off annoying warning
            
            %% input check (function may be called directly by user)
            
            if ~isscalar(BoxSize), error('BoxSize must be a scalar - to use an array, call VariableBoxSize_Analysis instead!!');end
            if BoxSize > obj.width || BoxSize > obj.height
                BoxSize = min(obj.width, obj.height);
                cprintf('*[1 .3 0]',['Impossible to run the analysis with the desired BoxSize, as it is bigger than at least one of the dimensions of the field of view. The analysis will be run on the maximum possible BoxSize, i.e. BoxSize = ',num2str(BoxSize)]);
            end
            if mod(BoxSize,2) %if odd BoxSize reduce of 1
                cprintf('*[1 .3 0]','BoxSize can''t be odd, reducing it of a unit.');
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
            
            N_Boxes_row = floor(obj.height / BoxSize);   %number of boxes that fit (vertically) in the field of view given BoxSize
            N_Boxes_col = floor(obj.width  / BoxSize);    %number of boxes that fit (horizontally) in the field of view given BoxSize
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
            
            col_span = N_Boxes_col * BoxSize;     %width of total DDM-analysed region (multiple of BoxSize)
            row_span = N_Boxes_row * BoxSize;    %height of total DDM-analysed region (multiple of BoxSize)
            
            %find where to place (in the horizontal direction) the [row_span, col_span] region to be DDM-analysed
            for i = obj.width - col_span + 1 : -1 : 1  %sneaky allocation
                dummy(i) = sum(obj.col_sum_std_fs(i:i+col_span-1));
            end
            [~,col_offset] = max(dummy);
            clear dummy
            %find where to place (in the vertical direction) the [row_span, col_span] region to be DDM-analysed
            for i = obj.height - row_span + 1 : -1 : 1  %sneaky allocation
                dummy(i) = sum(obj.row_sum_std_fs(i:i+row_span-1));
            end
            [~,row_offset] = max(dummy);
            clear dummy
            cprintf('*[0 .5 0]','Done!')
            
            %% select good boxes based on thresholding of standard deviation
            
            ind_good_boxes = find_good_boxes(obj.std_fs(row_offset:row_offset+row_span-1, col_offset:col_offset+col_span-1), BoxSize);
            
            %% saving in each Box the relative bit of std_fs
            
            fprintf('\n\t\tSaving portion of standard deviation of pixel intensity in time relative to each Box... ');
            for i = 1:N_Boxes_row
                for j = 1:N_Boxes_col
                    
                    ii = sub2ind([N_Boxes_row, N_Boxes_col], i, j);
                    obj.Results(iir).Box(ii).std_fs = obj.std_fs(row_offset + (i-1)*BoxSize : row_offset + (i)*BoxSize -1 ,...
                        col_offset + (j-1)*BoxSize : col_offset + (j)*BoxSize -1);
                end
            end
            cprintf('*[0 .5 0]','Done!')
            
            %% running DDM on each region of the video, storing Iqtaus in relative Box
            
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
            
            
            %% running fit_Iqtau_ourcilia on each Box.Iqtau
            
            fprintf('\n\t\tFitting Iqtau matrix for each Box... ');
            fprintf('\n\t\t\tFitting Box ');
            for ii = 1:N_Boxes
                
                fprintf('%.6d / %.6d',ii,N_Boxes);
                [obj.Results(iir).Box(ii).Frequency, obj.Results(iir).Box(ii).Damping, obj.Results(iir).Box(ii).Amplitude, obj.Results(iir).Box(ii).Offset, obj.Results(iir).Box(ii).GOF] = fit_Iqtau_ourcilia(obj.Results(iir).Box(ii).Iqtau, max_mode_fitted, obj.max_tau_fitted);
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
            end
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % deletes "fitting box"
            
            cprintf('*[0 .5 0]','Done!');
            
            
            %% saving parameters in results
            
            obj.Results(iir).row_offset = row_offset;
            obj.Results(iir).col_offset = col_offset;
            obj.Results(iir).ind_good_boxes = ind_good_boxes;
            
        end
        
        %% interactive creation of kymograph
        function obj = kymograph(obj,N_kymographs)
            
            if nargin==1,
                N_kymographs = 1;
            end
            
            global fs;
            if size(fs) ~= [obj.height, obj.width, obj.NumberOfFrames]
                cprintf('*[1 .3 0]','It looks like the video in the memory is not the same as specified by the instance of the class. Loading the proper video...');
                fs = load_movie(obj);
            end
            
            fprintf('\nAll previous kimographs for this video will be deleted!\n');
            
            obj.Kymograph = repmat(struct('ind',[], 'kymo',[]),[N_kymographs,1]);
            
            hf = figure(11001);
            imshow(fs(:,:,1),[]);
            
            for ki = 1:N_kymographs
                [col_kind, row_kind, ffprofile] = improfile;
                col_kind = round(col_kind);
                row_kind = round(row_kind);
                obj.Kymograph(ki).ind = sub2ind([obj.height, obj.width], row_kind, col_kind);
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
        function [obj, movie_correctly_read] = load_movie(obj)
            
            movie_correctly_read = false;   %control variable for failsafe movie reading
            
            global fs;
            if isprop(obj,'Filepath')
                fullpath = fullfile(obj.Filepath,obj.Filename);
            else
                fullpath = obj.Filename;
            end
            
            [~,~,ext] = fileparts(fullpath);
            
            if strcmp(ext,'.movie')
                try
                    mo = moviereader(fullpath);
                catch
                    cprintf('Red','The movie is not present at the path specified ');
                    return
                end
                
                obj.IsMovieCorrupted = false;
                
                while movie_correctly_read == false %try to read the movie until it works
                    try
                        fs = mo.read([obj.FirstFrame obj.LastFrame]);  %try to read the movie
                        movie_correctly_read = true;        %if it works put flag to 1
                        
                    catch err                               %if it doesn't work, then
                        if strcmp(err.identifier,'moviereader:WrongMagicAtFrame')
                            obj.LastFrame = str2double(err.message)-1;  %try to read again, ditching the last frame (I guess I could use a "finding" algorithm here)
                        else
                            obj.LastFrame = obj.LastFrame - 1;
                        end
                        obj.IsMovieCorrupted = true;        %set the flag that the movie was corrupted
                        
                    end                                     %and try again
                    
                end
                
                
                fs(:,:,end) = []; %ditches the last frame anyway, which may still be corrupted
                obj.LastFrame = obj.LastFrame - 1;
                
                if obj.IsMovieCorrupted
                    cprintf('Red','The movie is corrupted ')
                end
                
                
                
                obj.FrameRate   = mo.FrameRate;
                obj.NumberOfFrames = obj.LastFrame - obj.FirstFrame + 1;
                
                
                
            else
                fs = avi2greyscaleframestack(fullpath);
                fs = fs(:,:,obj.FirstFrame:obj.LastFrame);
            end
            
            if obj.NumberOfFrames ~= size(fs,3)
                error('mistake here')
            end
            
        end
        
        %% parse filename for magnification string
		% CAREFUL this has to be modified according to which camera you're using
        function obj = set_magnification(obj)
            
            if ~isempty(regexpi(obj.Filename,'4X'))
                if ~isempty(regexpi(obj.Filename,'1.5X'))
                    obj.Magnification = '6X';
                    obj.px2mum = 0.9727;
                else
                    obj.Magnification = '4X';
                    obj.px2mum = 1.459;
                end
            elseif ~isempty(regexpi(obj.Filename,'10X'))
                if ~isempty(regexpi(obj.Filename,'1.5X'))
                    obj.Magnification = '15X';
                    obj.px2mum = 0.3893;
                else
                    obj.Magnification = '10X';
                    obj.px2mum = 0.584;
                end
            elseif ~isempty(regexpi(obj.Filename,'20X'))
                if ~isempty(regexpi(obj.Filename,'1.5X'))
                    obj.Magnification = '30X';
                    obj.px2mum = 0.195;
                else
                    obj.Magnification = '20X';
                    obj.px2mum = 0.292;
                end
            elseif ~isempty(regexpi(obj.Filename,'30X'))
                if ~isempty(regexpi(obj.Filename,'1.5X'))
                    obj.Magnification = '45X';
                    obj.px2mum = 0.13;
                else
                    obj.Magnification = '30X';
                    obj.px2mum = 0.195;
                end
            elseif ~isempty(regexpi(obj.Filename,'40X'))
                if ~isempty(regexpi(obj.Filename,'1.5X'))
                    obj.Magnification = '60X';
                    obj.px2mum = 0.0973;
                else
                    obj.Magnification = '40X';
                    obj.px2mum = 0.146;
                end
            elseif ~isempty(regexpi(obj.Filename,'60X'))
                if ~isempty(regexpi(obj.Filename,'1.5X'))
                    obj.Magnification = '90X';
                    obj.px2mum = 0.0647;
                else
                    obj.Magnification = '60X';
                    obj.px2mum = 0.097;
                end
            else
                cprintf('*[1 .3 0]','\n Couldn''t find magnification string in namefile.');
            end
            
        end
        
        %% function that gathers the results (does means etc)
        function obj = gather_results(obj)
            
            for ii = 1:numel(obj.Results) %create a vector with median frequency of each box: as long as the number of boxes for each boxsize
                
                if ~isfield(obj.Results(ii),'ind_good_boxes') || (isfield(obj.Results(ii),'ind_good_boxes') && isempty(obj.Results(ii).ind_good_boxes) )
                    row_span = floor(obj.height / obj.Results(ii).BoxSize) * obj.Results(ii).BoxSize;   %number of boxes that fit (vertically) in the field of view given BoxSize
                    col_span = floor(obj.width  / obj.Results(ii).BoxSize) * obj.Results(ii).BoxSize;
                    ind_good_boxes = find_good_boxes(obj.std_fs(obj.Results(ii).row_offset:obj.Results(ii).row_offset+row_span-1,obj.Results(ii).col_offset:obj.Results(ii).col_offset+col_span-1), obj.Results(ii).BoxSize);
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
                
                obj.Results(ii).MedianFrequencyVec = nanmedian(temp_Frequency_mat); % bc = box_counter %% IN HERTZ NOW!!!!
                obj.Results(ii).AverageDamping_vs_q = trimmean( temp_Damping_mat ,20 ,2);   %% IN HERTZ NOW!!!!
                obj.Results(ii).AverageAmplitude_vs_q = trimmean( temp_Amplitude_mat, 10 ,2);
                
                obj.Results(ii).MedianDamping_vs_q = nanmedian( temp_Damping_mat ,2);   %% IN HERTZ NOW!!!!
                obj.Results(ii).MedianAmplitude_vs_q = nanmedian( temp_Amplitude_mat ,2);
                
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
                bsz = mean(floor([obj.height, obj.width]./size(obj.SAVAlike.frequency_map)));
                obj.SAVAlike.ind_good_bins = find_good_boxes(obj.std_fs,bsz);
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
            
            
            rect = [1 1 obj.width obj.height];
            
            if flag_ROI %select ROI
                figure(11002);
                imagesc(obj.std_fs); colormap jet;
                
                rect = round(getrect(gcf));
                
                rect(3) = bsz*ceil(rect(3)/bsz);
                rect(4) = bsz*ceil(rect(4)/bsz);
                
                fs = fs(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1,:);
            end
            
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
            
            if nargin < 2
                channel = 0; % select TC1
            end
            channel = channel + 2; %to take into account the time vector column and C notation
            
            TemperatureFilename = [obj.Filename(1:end-5),'temperature'];
            if isprop(obj,'Filepath')
                TemperatureFilename = fullfile(obj.Filepath,TemperatureFilename);
            end
            if isempty(dir(TemperatureFilename))
                fprintf('\n');
                cprintf('*[1 .3 0]',regexprep(['Temperature file not found at ',TemperatureFilename,'.'],'\','\\\\'));
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
        
        
    end
    
    
end

%% select good boxes based on thresholding of standard deviation

function mask = find_good_boxes(IM, bsz)

IM = mat2gray(IM);      %converts to [0 1]
[height, width] = size(IM);
binningmap = col2im(repmat(1:floor(width*height/bsz^2),bsz*bsz,1), bsz.*[1 1], bsz*[floor(height/bsz) floor(width/bsz)],'distinct');
bIM = reshape( accumarray(binningmap(:),IM(:)), height/bsz, width/bsz) ./ bsz.^2;
abIM = imadjust(bIM);
mask = im2bw(abIM,min(1,graythresh(bIM)));

if bsz < 128
    mask = bwmorph(mask,'clean'); % maybe want to think about this again
end

end
