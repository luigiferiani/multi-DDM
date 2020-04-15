% Thanks to Dr Elvis Pandzic from the University of New South Wales, 
% for sharing his code to read bioformat videos.
% For the moment only single-channel videos are supported, but this will be
% improved over time

classdef OmeReaderWrapper < matlab.mixin.Copyable
    
    properties
        Filename            %path of the video file
        FirstFrame          %first frame to analyse
        LastFrame           %last frame to analyse
        Height              %rows of video, (px)
        Width               %cols of video, (px)
        TotalNumberOfFrames %N_frames as written in the movie header
        NumberOfFrames      %number of frames to analyse (LastFrame-FirstFrame+1)
        FrameRate           %
        Magnification       %for px2mum conversions
        px2mum              %pixels to micron conversion
    end % properties
    
    properties (SetAccess = private, Hidden = true)
        vo
        specs
    end
    
    methods
        % constructor
        function obj = OmeReaderWrapper(filename, firstframe, lastframe)
            
            obj.Filename = filename;
            obj.vo = bfGetReader(filename);
            
            % read metadata
            obj.read_metadata;
            
            % number of frames
            
            if nargin >= 3 && ~isempty(lastframe)
                obj.LastFrame = min(obj.TotalNumberOfFrames, lastframe);
            else
                obj.LastFrame = obj.TotalNumberOfFrames;
            end % if
            
            if nargin >= 2 && ~isempty(firstframe)
                obj.FirstFrame = max(1, firstframe);
            else
                obj.FirstFrame = 1;
            end % if
            
            obj.NumberOfFrames = obj.LastFrame - obj.FirstFrame + 1;
            
        end % function
        
        % metadata
        function obj = read_metadata(obj)
            %read_metadata this uses the code by E. Pandzic, then
            %translates it for the wrapper
            
            % read OME metadata
            omeMeta = obj.vo.getMetadataStore();
            
            % parse px->um
            try
                obj.specs.voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value; % in µm
                obj.specs.voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value; % in µm
                obj.specs.pixelSize = (str2double(obj.specs.voxelSizeX)+str2double(obj.specs.voxelSizeY))/2;
                obj.specs.is_RWU_space = true;
            catch
                obj.specs.is_RWU_space = false;
                obj.specs.pixelSize = 1;
            end
            
            % parse frame -> s
            obj.specs.sizex = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels;
            obj.specs.sizey = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels;
            obj.specs.TotImages = omeMeta.getPixelsSizeT(0).getValue(); %number of images
            try
                obj.specs.timeFrame = double(omeMeta.getPlaneDeltaT(0,1).value)-double(omeMeta.getPlaneDeltaT(0,0).value);
                obj.specs.is_RWU_time = true;
            catch
                obj.specs.is_RWU_time = false;
                obj.specs.timeFrame = 1;
            end
            
            obj.specs.numChannels = omeMeta.getChannelCount(0);
            obj.specs.stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue();
            %specs.NA = omeMeta.getObjectiveLensNA(0,0).doubleValue();
            %specs.laser = omeMeta.getChannelExcitationWavelength(0, 0).value().doubleValue();
            
            % translate now for wrapper
            obj.Width = obj.specs.sizex;
            obj.Height = obj.specs.sizey;
            obj.TotalNumberOfFrames = obj.specs.TotImages;
            if obj.specs.is_RWU_time
                obj.FrameRate = 1/obj.specs.timeFrame;
            end % else empty, for compatibility
            if obj.specs.is_RWU_space
                obj.px2mum = obj.specs.pixelSize;
            end % else empty, for compatibility
            
        end % function
        
        % reader
        function fs = read(obj, frame_boundaries)
            %read frameboundaries allows to bypass the firstframe lastframe
            %inputs in the class constructor. only used if you want to read
            %a sub-subset of frames, e.g. the 1st second in DDM_Analysis
            %for motion estimation
            if nargin < 2 || isempty(frame_boundaries) 
                frame_boundaries = [obj.FirstFrame, obj.LastFrame];
            end %if
            number_of_frames = diff(frame_boundaries) + 1;
            
            ff = obj.vo.bfGetPlane(1); % read first frame to get type
            fs = zeros(obj.Height, obj.Width, number_of_frames, class(ff));
            
            if obj.specs.numChannels == 1
                
                for dat = frame_boundaries(1) : frame_boundaries(2)
                    fs(:,:,dat) = obj.vo.bfGetPlane(dat);
                end % for
            else
                error('Only one-channel videos are supported at the moment')
            end % if
            
%             if obj.specs.numChannels == 1
%                 
%                 for dat =1:obj.TotalNumberOfFrames
%                     plane = bfGetPlane(reader,dat);
%                     image_data(:,:,dat) = single(plane);
%                 end
%             elseif obj.specs.numChannels == 2
%                 for dat=1:2*obj.TotalNumberOfFrames
%                     plane=bfGetPlane(reader,dat);
%                     data(:,:,dat)=single(plane);
%                 end
%                 
%                 for j=1:round(size(data,3)/2)
%                     image_data(:,:,1,j)=data(:,:,1+(j-1)*2);
%                     image_data(:,:,2,j)=data(:,:,j*2);
%                 end
%                 
%             elseif obj.specs.numChannels == 3
%                 for dat=1:3*obj.TotalNumberOfFrames
%                     plane=bfGetPlane(reader,dat);
%                     data(:,:,dat)=single(plane);
%                 end
%                 
%                 for j=1:round(size(data,3)/3)
%                     image_data(:,:,1,j)=data(:,:,1+(j-1)*3);
%                     image_data(:,:,2,j)=data(:,:,2+(j-1)*3);
%                     image_data(:,:,3,j)=data(:,:,j*3);
%                 end
%             end
            
        end % function
        
    end % methods
    
end % classdef

