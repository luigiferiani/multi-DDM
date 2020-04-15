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
        vid
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
        function fs = read(obj)
            user=memory;
            user.MemUsedMATLAB/(1024)^2
            % simply call the moviereader method
            fs = obj.vo.read([obj.FirstFrame, obj.LastFrame]);
            user=memory;
            user.MemUsedMATLAB/(1024)^2
        end % function
        
    end % methods
    
end % classdef

