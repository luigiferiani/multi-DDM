function vid = select_video_reader(filename, varargin)

is_moviereader = any(endsWith(filename, ".movie"));
is_videoreader = any(endsWith(filename, [".avi", ".mp4"]));
is_omereader   = any(endsWith(filename, [".czi", ".lif", ".oif"]));


if is_moviereader
    
    vid = MovieReaderWrapper(filename, varargin{:});
    
elseif is_videoreader
    
    vid = VideoReaderWrapper(filename, varargin{:});
    
elseif is_omereader
    
    vid = OmeReaderWrapper(filename, varargin{:});

else
    
    % duck typing
    try % videoreader
        vid = VideoReaderWrapper(filename, varargin{:});
    catch
        try % omereader
            vid = OmeReaderWrapper(filename, varargin{:});
        catch
            error("Movie format not supported at this stage")
        end % try omereader
    end % try videoreader
    
end % if


end % function