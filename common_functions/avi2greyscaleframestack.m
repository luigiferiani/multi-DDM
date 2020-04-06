function [ frame_stack, video ] = avi2greyscaleframestack( filename )
%avi2framestack Converts avi files to a 3D matrix


video = VideoReader(filename);
video_data = video.read;
frame_stack = zeros(video.Height, video.Width, video.NumberOfFrames,'uint8');
VideoFormat = video.VideoFormat;

switch VideoFormat
    case 'RGB24'
        for i=1:video.NumberOfFrames
            frame_stack(:,:,i) = rgb2gray(video_data(:,:,:,i));
        end
        
    case 'Grayscale'
        frame_stack = squeeze(video_data);
        
    otherwise
        warning('Unexpected VideoFormat. Modify the function accordingly.');
end

end

