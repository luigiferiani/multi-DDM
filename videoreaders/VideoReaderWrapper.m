classdef VideoReaderWrapper < matlab.mixin.Copyable
    
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
        function obj = VideoReaderWrapper(filename, firstframe, lastframe)
            
            obj.Filename = filename;
            obj.vo = VideoReader(filename);
            
            % read metadata
            for foo = {'Height', 'Width', 'FrameRate'}
                prop = foo{1};  % tostring
                obj.(prop) = obj.vo.(prop);
            end % for
            % number of frames
            obj.TotalNumberOfFrames = obj.vo.Duration * obj.FrameRate;
            
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
%             user=memory;
%             user.MemUsedMATLAB/(1024)^2
            % reset video if necessary
            if obj.vo.CurrentTime > 0
                obj.vo = VideoReader(obj.Filename);
            end % if
            % how to treat each frame
            switch obj.vo.VideoFormat
                case 'RGB24'
                    squeezing_func = @rgb2gray;
                case 'Grayscale'
                    squeezing_func = @squeeze;
                otherwise
                    warning('Unexpected VideoFormat. Modify the function accordingly.');
            end % switch
            
            % shift video beginning if needed
            if obj.FirstFrame > 1
                obj.vo.CurrentTime = (obj.FirstFrame - 1) / obj.FrameRate;
            end % if
            
            % read
            %             global fs
            fs = zeros(obj.Height, obj.Width, obj.NumberOfFrames, 'uint8');
            for fc = 1:obj.NumberOfFrames
                if obj.vo.hasFrame
                    fs(:,:,fc) = squeezing_func(obj.vo.readFrame);
                else
                    error('out of Frames, off by one error?')
                end
            end % for
%             user=memory;
%             user.MemUsedMATLAB/(1024)^2
        end % function
        
    end % methods
    
end % classdef

