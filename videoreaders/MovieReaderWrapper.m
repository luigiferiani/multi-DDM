classdef MovieReaderWrapper < matlab.mixin.Copyable
    
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
    end
    
    methods
        % constructor
        function obj = MovieReaderWrapper(filename, firstframe, lastframe)
            
            obj.Filename = filename;
            obj.vo = moviereader(filename);
            
            % metadata
            obj.Height = double(obj.vo.height);
            obj.Width  = double(obj.vo.width);            
            
            % number of frames
            obj.TotalNumberOfFrames = obj.vo.NumberOfFrames;
            obj.FrameRate = obj.vo.FrameRate;
            
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
        
        % reader
        function fs = read(obj, frame_boundaries)
            %read frameboundaries allows to bypass the firstframe lastframe
            %inputs in the class constructor. only used if you want to read
            %a sub-subset of frames, e.g. the 1st second in DDM_Analysis
            %for motion estimation
            if nargin < 2 || isempty(frame_boundaries) 
                frame_boundaries = [obj.FirstFrame, obj.LastFrame];
            end %if

            % simply call the moviereader method
            fs = obj.vo.read(frame_boundaries);

        end % function
        
    end % methods
    
end % classdef

