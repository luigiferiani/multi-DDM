% A very Basic Matlab interface to the .movie binary files
classdef moviereader < handle
    
    % movie properties and data
    properties
        Directory % The subdirectory that contains the video
        Filename  % The filename of the video
        NumberOfFrames % The number of frames contained in the video
        FrameRate % The number of frames recorded in a second
        width     % The image width, in pixels
        height    % The image height, in pixels
        
        
        
    end
    
    properties (SetAccess=private, Hidden)
        offset_header; % bytes before the data
        version;
        camera_type;
        endian;
        length_data;
        length_header;
        data_depth;
        image_bytes;
        total_bytes;
        length_in_words;
        first_frame_timestamp_sec;
        magic;
    end
    
    % functions to read and display the movie
    methods
        function movie=moviereader(moviename,varargin)
            % Function that loads the video and initialize its
            % structure, extracting the infos stored into the header of the
            % movie file.
            %
            % USAGE: movie=moviereader(filename);
            %        movie=moviereader(filename,directory);
            
            datadir='';
            if nargin>1
                datadir=varargin{1};
                if ~strcmp(datadir(end),'/');
                    datadir=[datadir,'/'];
                end
            end
            
            
            movie.Directory=datadir;
            movie.Filename=moviename;
            
            fid=fopen([datadir,moviename],'r');
            % Camera Constants
            % % Finding begining of the frame header need to store the offset
            % % The ASCII part is only stored once, before 1st frame's header
            offset=0;
            word=0;
            
            CAMERA_MOVIE_MAGIC = hex2dec('496d6554'); %'TemI'
            CAMERA_TYPE_IIDC = 1;
            CAMERA_TYPE_ANDOR= 2;
            CAMERA_TYPE_XIMEA= 3;
            
            % % find offset of the first header  store in offset_header
            while ( word ~= CAMERA_MOVIE_MAGIC )
                
                % % is there another word in the file?
                % add is not in a MB stop !!!
                status = fseek( fid, offset, -1 ); % from begining (-1)
                if ( status == -1)
                    fprintf('No Magic word found in all of the file \n');
                    return
                end
                word = fread( fid, 1,'*uint32' ) ;
                offset = offset + 1  ;
                if (offset >= 1000000)
                    fprintf('The potential ASCII header is too long. \n No "TemI" found, bailing out. \n Old data type( pre January 2012) ? \n');
                    return
                end
                % %  If it is magic word note the offset
                if ( word == CAMERA_MOVIE_MAGIC )
                    offset = offset -1 ; % point back to magic
                    fseek( fid, offset, -1 ) ;
                    break;
                end
            end
            
            % This variable is important make sure to keep it for as long as you
            % read from this file - additional offset before first header
            movie.offset_header = offset;
            
            
            % Common stuff for all camera types
            common_info = fread(fid,6,'*uint32'); %read 24 bytes
            movie.version = common_info(2);
            movie.camera_type = common_info(3);
            
            if (movie.version == 1)
                movie.endian=common_info(4); % BigEndian,16bit==2, LittleEndian,16bit==1, 8bit==1, ??12bit==4??
                movie.length_header = common_info(5);
                movie.length_data = common_info(6);
            elseif (movie.version == 0)
                %need to hardwire BigEndian/LittleEndian
                movie.endian = 2; % BigEndian==2
                movie.length_header = common_info(4) ;
                movie.length_data = common_info(5)  ;
                fseek(fid,-4,'cof') ; % go back one extra word read into the common_info
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
            % Camera specific settings
            switch movie.camera_type
                
                %%%%%% IIDC %%%%%%%%
                case CAMERA_TYPE_IIDC
                    %                     disp('IIDC (Grasshopper) camera')
                    fseek( fid, 24, 0 ); % skip next 24 bytes to c_timestamp 44 position
                    %                     if (movie.endian == 2)
                    %                         c_timestamp([2 1],1) = fread(fid,2,'*uint32'); %shouldn't this be just a 'uint64'?
                    %                     else
                    %                         c_timestamp([1 2],1) = fread(fid,2,'*uint32');
                    %                     end
                    movie.first_frame_timestamp_sec = double(fread(fid,1,'*uint64'))/10^6;
                    if (movie.version == 0)
                        data_shape([1 2 3 4 5 6 7 8]) = fread(fid,8,'*uint32');
                        data_shape(9) = fread(fid,1,'*uint64'); %total data stored
                    elseif (movie.version == 1)
                        data_shape([1 2 3 4 5 6 7 8 9 10]) = fread(fid,10,'*uint32');
                        data_shape(11) = fread(fid,1,'*uint64'); %total data i_total_bytes stored
                    end
                    % columns are:  (v.1 only i_size_x_max; i_size_y_maxc_size_x;) c_size_y; c_pos_x; c_pos_y; c_pixnum; c_stride;
                    % c_data_depth; c_image_bytes; c_total_bytes
                    
                    % go to the c_timestamp in the next header :
                    if (movie.version == 0)
                        % (44 bytes from begining of that header to position before c_timestamp):
                        timestamp_pos = 44;
                    elseif (movie.version == 1)
                        % (48 bytes from begining of that header to position befor c_timestamp):
                        timestamp_pos = 48;
                    end
                    
                    fseek( fid, movie.offset_header + movie.length_header + movie.length_data +timestamp_pos, -1 );
                    %                     if (movie.endian == 2)
                    %                         c_timestamp([2 1],2) = fread(fid,2,'*uint32');
                    %                     else
                    %                         c_timestamp([1 2],2) = fread(fid,2,'*uint32');
                    %                     end
                    
                    
                    %                     movie.FrameRate = 1/( ( double(c_timestamp(1,2)) - double(c_timestamp(1,1))) /10^6 );  %find out frame rate (f/s)
                    movie.FrameRate = 1/(double(fread(fid,1,'*uint64'))/10^6 - movie.first_frame_timestamp_sec);
                    movie.width = data_shape(3);
                    movie.height= data_shape(4);
                    
                    
                    
                    movie.data_depth = data_shape(9);  %this is number of bytes 8 or 16
                    
                    movie.image_bytes = data_shape(10);
                    movie.total_bytes = data_shape(11);
                    
                    if (movie.data_depth == 16)
                        movie.length_in_words = floor(movie.image_bytes/2); %each word is 2 bytes;
                    elseif (movie.data_depth == 8)
                        movie.length_in_words = movie.image_bytes;
                    elseif (movie.data_depth == 12)
                        movie.length_in_words = floor(movie.image_bytes/1.5); %each word is an abomination of 12 bits, aka 1.5 bytes;
                    else
                        fprintf('Bit depth neither 8 nor 12 nor 16 - not programmed for this.\n');
                        return;
                    end
                    
                    %%%%%% ANDOR %%%%%%%%
                case CAMERA_TYPE_ANDOR
                    movie.endian = 1; % LittleEndian==1
                    disp('ANDOR camera')
                    
                    c_timestamp(1,[1,2]) = fread(fid,2,'*uint64');       % time of the first image
                    movie.first_frame_timestamp_sec = double(c_timestamp(1,1)) + double(c_timestamp(1,2))/10^9;
                    
                    if (movie.version == 1)
                        data_shape([1 2 3 4 5 6 7 8]) = fread(fid,8,'*uint32');
                        timestamp_pos = 24;
                    end
                    
                    fseek( fid, movie.offset_header + movie.length_header + movie.length_data +timestamp_pos, -1 );
                    c_timestamp(2,[1,2]) = fread(fid,2,'*uint64');
                    
                    movie.FrameRate = 1/( ( double(c_timestamp(2,1)) + double(c_timestamp(2,2))/10^9) - ( double(c_timestamp(1,1)) + double(c_timestamp(1,2))/10^9));  %find out frame rate (f/s)
                    
                    movie.width = data_shape(4) - data_shape(3) +1;
                    movie.height= data_shape(6) - data_shape(5) +1;
                    
                    movie.data_depth = 8*movie.length_data/movie.height/movie.width;
                    
                    
                    %%%%%% XIMEA %%%%%%%%
                case CAMERA_TYPE_XIMEA
                    disp('XIMEA camera')
                    fread(fid,100,'*char'); %camera name... I don't need that, I just skip
                    dim_char=100;
                    fseek( fid, 4, 0 ); % skip next 4 bytes to go to timestamp
                    
                    c_timestamp(1,[1,2]) = fread(fid,2,'*uint64');       % time of the first image
                    movie.first_frame_timestamp_sec = double(c_timestamp(1,1)) + double(c_timestamp(1,2))/10^9;
                    
                    if (movie.version == 1)
                        data_shape([1 2 3 4 5 6]) = fread(fid,6,'*uint32');
                        timestamp_pos = 28+dim_char;
                    end
                    
                    fseek( fid, movie.offset_header + movie.length_header + movie.length_data +timestamp_pos, -1 );
                    c_timestamp(2,[1,2]) = fread(fid,2,'*uint64');
                    
                    
                    movie.FrameRate = 1/( ( double(c_timestamp(2,1)) + double(c_timestamp(2,2))/10^9) - ( double(c_timestamp(1,1)) + double(c_timestamp(1,2))/10^9));  %find out frame rate (f/s)
                    movie.width = data_shape(3);
                    movie.height= data_shape(4);
                    
                    movie.length_data;
                    movie.data_depth = 8*movie.length_data/movie.height/movie.width;  %this is number of bytes 8 or 16
                    
                    %checking if the movie contains full or cropped frames
                    if (data_shape(3) == movie.width) && (data_shape(4) == movie.height)
                        movie.image_bytes = movie.width*movie.height;
                    else
                        movie.image_bytes = data_shape(3)*data_shape(4);            % this is if the image has been cropped
                    end
                    
                    if (movie.data_depth == 16)
                        movie.length_in_words = floor(movie.image_bytes/2); %each word is 2 bytes;
                    elseif (movie.data_depth == 8)
                        movie.length_in_words = movie.image_bytes;
                    else
                        fprintf('data neither 8 not 16 bit - not programmed for this.\n');
                        return;
                    end
                    
            end
            fseek(fid, 0, 'eof');
            file_size = ftell(fid) ;
            
            %position at the begining of the binary header+frame data
            fseek( fid, movie.offset_header , -1 );
            
            nn = (file_size - movie.offset_header) / (movie.length_header + movie.length_data);
            if (floor(nn) ~= nn)
                fprintf('something is wrong in the calculation nn needs to be an integer \n');
            end
            
            movie.NumberOfFrames=double(nn);
            fclose(fid);
            
            movie.magic = CAMERA_MOVIE_MAGIC;
            
        end
        
        %                         function  [data]=read(movie,varargin)
        function [IM, timestamp, header]=read(movie,varargin)
            % Function that loads frames contained in the movie
            % movie is a moviereader object, previously created
            %
            % USAGE: frames=movie.read;              -load all the frames
            %        frames=movie.read(num);         -load the nth frame
            %        frames=movie.read([num1,num2]); -load the frames num1:num2
            %
            % ALT. USAGE: frames=read(movie);
            %             frames=read(movie,num);
            %             frames=read(movie,[num1,num2]);
            %
            %optional input parameter: frames to be loaded
            if nargin==2
                range=varargin{1};
                if length(range)>2;
                    disp('Cannot load non contiguous frames. Call load iteratively');
                    return
                end
            elseif nargin==1;
                range=[1,movie.NumberOfFrames];
            end
            
            if length(range)==2
                if range(2)>movie.NumberOfFrames
                    disp(['The video contains ',num2str(movie.NumberOfFrames),' frames only'])
                    range(2)=movie.NumberOfFrames;
                end
            elseif length(range)==1
                if range(1)>movie.NumberOfFrames
                    disp(['The video contains ',num2str(movie.NumberOfFrames),' frames only'])
                    range(1)=movie.NumberOfFrames;
                end
            end
            
            % constants
            CAMERA_TYPE_IIDC = 1;
            CAMERA_TYPE_ANDOR= 2;
            CAMERA_TYPE_XIMEA= 3;
            
            % bytes to skip %% THIS LINE HAS TO BE MODIFIED IF THE FRAME
            % SIZE IS NOT CONSTANT sth like while counter ~=
            % frame_I_want_to_read => read the header of each frame and
            % skip the image data
            offset=uint64(movie.offset_header)+ uint64(range(1)-1) * (uint64(movie.length_header) + uint64(movie.length_data));
            offset = double(offset); %to prevent errors on some computers
            
            if length(range)>1
                N_frames_to_load=diff(range)+1;
            else
                N_frames_to_load=1;
            end
            
            
            %             fmap=memmapfile([movie.Directory,movie.Filename],...
            %                 'Offset', offset,...
            %                 'Format', {'uint8' double([1, movie.length_header]) 'info';
            %                            ['uint',int2str(movie.data_depth)] double([movie.width, movie.height]) 'IM'},...
            %                 'Repeat', N_frames_to_load);
            switch movie.camera_type
                
                %%%%%% IIDC %%%%%
                case CAMERA_TYPE_IIDC
                    
                    if movie.data_depth == 12
                        
                        conv = repmat(uint16([2^8 2^4 1  1 2^8 2^4 ]'),[double(movie.width)*double(movie.height)/2, 1]); %conversion matrix
                        pxllist = repmat(1:double(movie.width)*double(movie.height),[3,1]);
                        pxllist = pxllist(:);
                        pxllist = pxllist + repmat([0 0 +1 -1 0 0 ],[1 double(movie.width)*double(movie.height)/2])';
                        
                        fid=fopen([movie.Directory,movie.Filename],'r');
                        fseek( fid, offset , -1 );
                        for i=1:N_frames_to_load
                            data(i).magic           = fread(fid,1,'*uint32'); %40
                            data(i).version         = fread(fid,1,'*uint32'); %44
                            data(i).type            = fread(fid,1,'*uint32'); %48
                            data(i).pixelmode       = fread(fid,1,'*uint32'); %52
                            data(i).length_header   = fread(fid,1,'*uint32'); %56
                            data(i).length_data     = fread(fid,1,'*uint32'); %60
                            data(i).i_guid          = fread(fid,1,'*uint64'); %64
                            data(i).i_vendor_id     = fread(fid,1,'*uint32'); %72
                            data(i).i_model_id      = fread(fid,1,'*uint32'); %76
                            data(i).i_video_mode    = fread(fid,1,'*uint32'); %80
                            data(i).i_color_coding  = fread(fid,1,'*uint32'); %84
                            data(i).i_timestamp_us  = fread(fid,1,'*uint64'); %88
                            data(i).i_size_x_max    = fread(fid,1,'*uint32'); %96
                            data(i).i_size_y_max    = fread(fid,1,'*uint32'); %100
                            data(i).i_size_x        = fread(fid,1,'*uint32'); %104
                            data(i).i_size_y        = fread(fid,1,'*uint32'); %108
                            data(i).i_pos_x         = fread(fid,1,'*uint32'); %112
                            data(i).i_pos_y         = fread(fid,1,'*uint32'); %116
                            data(i).i_pixnum        = fread(fid,1,'*uint32'); %120
                            data(i).i_stride        = fread(fid,1,'*uint32'); %124
                            data(i).i_data_depth    = fread(fid,1,'*uint32'); %128
                            data(i).i_image_bytes   = fread(fid,1,'*uint32'); %
                            data(i).i_total_bytes   = fread(fid,1,'*uint64');
                            data(i).i_brightness_mode = fread(fid,1,'*uint32'); %should try and solve the "union" thing here
                            data(i).i_brightness    = fread(fid,1,'*uint32');
                            data(i).i_exposure_mode = fread(fid,1,'*uint32');
                            data(i).i_exposure      = fread(fid,1,'*uint32');
                            data(i).i_gamma_mode    = fread(fid,1,'*uint32');
                            data(i).i_gamma         = fread(fid,1,'*uint32');
                            data(i).i_shutter_mode  = fread(fid,1,'*uint32');
                            data(i).i_shutter       = fread(fid,1,'*uint32');
                            data(i).i_gain_mode     = fread(fid,1,'*uint32');
                            data(i).i_gain          = fread(fid,1,'*uint32');
                            data(i).i_temperature_mode = fread(fid,1,'*uint32');
                            data(i).i_temperature   = fread(fid,1,'*uint32');
                            data(i).i_trigger_delay_mode = fread(fid,1,'*uint32');
                            data(i).i_trigger_delay = fread(fid,1,'*uint32');
                            data(i).i_trigger_mode  = fread(fid,1,'*int32');
                            data(i).i_avt_channel_balance_mode = fread(fid,1,'*uint32');
                            data(i).i_avt_channel_balance = fread(fid,1,'*int32');
                            
                            temp = fread(fid, movie.image_bytes*2, 'ubit4=>uint16' ,'b'); %reading 4bits words
                            % data(i).IM = reshape( uint16(accumarray(pxllist,temp.*conv)), [double(movie.width) double(movie.height)] );
                            data(i).IM = uint16(accumarray(pxllist,temp.*conv));
                            % sum of 4bitlong words belonging to the same pixel, appropriately multiplied by their "place"
                            
                            data(i).dummy = fread(fid, double(movie.total_bytes)-double(movie.image_bytes), '*uint32');
                            
                            if (i==1 && N_frames_to_load ~= 1) %allocation of all the subsequent frames
                                disp('allocating...');
                                % data(N_frames_to_load) = data(1);
                                data = repmat(data(1),[N_frames_to_load,1]);
                                disp('done');
                            end
                        end
                        
                        fclose(fid);
                        
                        for i=1:N_frames_to_load
                            if data(i).magic ~= movie.magic,
                                ME = MException('moviereader:WrongMagicAtFrame',num2str(i));
                                throw(ME);
                            end
                            data(i).timestamp_sec = double(data(i).i_timestamp_us)/10^6;
                            data(i).relative_timestamp_sec = data(i).timestamp_sec - movie.first_frame_timestamp_sec;
                        end
                        
                    else %use old version. should be faster
                        fmap=memmapfile([movie.Directory,movie.Filename],...
                            'Offset', offset,...
                            'Format', { 'uint32' 1 'magic';
                            'uint32' 1 'version';
                            'uint32' 1 'type';
                            'uint32' 1 'pixelmode';
                            'uint32' 1 'length_header';
                            'uint32' 1 'length_data';
                            'uint64' 1 'i_guid'; %this won't be readable, but I think 'chars' on MatLab are 2 bytes instead of 1 and I don't want to take any chance
                            'uint32' 1 'i_vendor_id';
                            'uint32' 1 'i_model_id';
                            'uint32' 1 'i_video_mode';
                            'uint32' 1 'i_color_coding';
                            'uint64' 1 'i_timestamp_us'; %in microseconds
                            'uint32' 1 'i_size_x_max';
                            'uint32' 1 'i_size_y_max';
                            'uint32' 1 'i_size_x';
                            'uint32' 1 'i_size_y';
                            'uint32' 1 'i_pos_x';
                            'uint32' 1 'i_pos_y';
                            'uint32' 1 'i_pixnum';
                            'uint32' 1 'i_stride';
                            'uint32' 1 'i_data_depth';
                            'uint32' 1 'i_image_bytes';
                            'uint64' 1 'i_total_bytes';
                            'uint32' 1 'i_brightness_mode';
                            'uint32' 1 'i_brightness';  %this could be a 'float' (in movie2tiff a 'union' was used, but I don't know how to handle that in MatLab).
                            'uint32' 1 'i_exposure_mode';
                            'uint32' 1 'i_exposure';  %this could be a 'float' (in movie2tiff a 'union' was used, but I don't know how to handle that in MatLab).
                            'uint32' 1 'i_gamma_mode';
                            'uint32' 1 'i_gamma';  %this could be a 'float' (in movie2tiff a 'union' was used, but I don't know how to handle that in MatLab).
                            'uint32' 1 'i_shutter_mode';
                            'uint32' 1 'i_shutter';  %this could be a 'float' (in movie2tiff a 'union' was used, but I don't know how to handle that in MatLab).
                            'uint32' 1 'i_gain_mode';
                            'uint32' 1 'i_gain';  %this could be a 'float' (in movie2tiff a 'union' was used, but I don't know how to handle that in MatLab).
                            'uint32' 1 'i_temperature_mode';
                            'uint32' 1 'i_temperature';  %this could be a 'float' (in movie2tiff a 'union' was used, but I don't know how to handle that in MatLab).
                            'uint32' 1 'i_trigger_delay_mode';
                            'uint32' 1 'i_trigger_delay';  %this could be a 'float' (in movie2tiff a 'union' was used, but I don't know how to handle that in MatLab).
                            'int32'  1 'i_trigger_mode';
                            'uint32' 1 'i_avt_channel_balance_mode';
                            'int32'  1 'i_avt_channel_balance';
                            
                            %                         ['uint',int2str(movie.data_depth)] double(movie.total_bytes) 'IM'},...
                            ['uint',int2str(8*ceil(movie.data_depth/8))] [double(movie.width) double(movie.height)] 'IM';
                            %                         ['uint',int2str(movie.data_depth)] double(movie.total_bytes)-double(movie.image_bytes) 'dummy'},...
                            ['uint',int2str(8*ceil(movie.data_depth/8))] (double(movie.total_bytes)-double(movie.image_bytes))/(double(movie.data_depth)/8) 'dummy'},...
                            'Repeat', N_frames_to_load);
                        
                        data=fmap.data;
                        clear fmap
                        
                        
                        for i=1:N_frames_to_load
                            if data(i).magic ~= movie.magic,
                                ME = MException('moviereader:WrongMagicAtFrame',num2str(i));
                                throw(ME);
                            end
                            data(i).timestamp_sec = double(data(i).i_timestamp_us)/10^6;
                            data(i).relative_timestamp_sec = data(i).timestamp_sec - movie.first_frame_timestamp_sec;
                            data(i).IM = swapbytes(data(i).IM);
                        end
                    end
                    
                    %%%%% ANDOR %%%%%
                case CAMERA_TYPE_ANDOR
                    fmap=memmapfile([movie.Directory,movie.Filename],...
                        'Offset', offset,...
                        'Format', { 'uint32' 1 'magic';
                        'uint32' 1 'version';
                        'uint32' 1 'type'; % Camera type
                        'uint32' 1 'pixelmode'; % Pixel mode
                        'uint32' 1 'length_header'; % Header data in bytes ( Everything except image data )
                        'uint32' 1 'length_data'; % Total data length in bytes;
                        'uint64' 1 'a_timestamp_sec';
                        'uint64' 1 'a_timestamp_nsec';
                        'int32'  1 'a_x_size_max'; % Sensor size
                        'int32'  1 'a_y_size_max';
                        'int32'  1 'a_x_start'; % Selected size and positions
                        'int32'  1 'a_x_end';
                        'int32'  1 'a_y_start';
                        'int32'  1 'a_y_end';
                        'int32'  1 'a_x_bin';
                        'int32'  1 'a_y_bin';
                        'int32'  1 'a_ad_channel'; % ADC
                        'int32'  1 'a_amplifier'; % EM or classical preamplifier
                        'int32'  1 'a_preamp_gain'; % Preamplifier gain
                        'int32'  1 'a_em_gain'; % EM gain
                        'int32'  1 'a_hs_speed'; % HS speed
                        'int32'  1 'a_vs_speed'; % VS speed
                        'int32'  1 'a_vs_amplitude'; % VS amplitude
                        'single' 1 'a_exposure'; % Exposure time in seconds
                        'int32'  1 'a_shutter'; % Shutter
                        'int32'  1 'a_trigger'; % Trigger
                        'int32'  1 'a_temperature'; % Temperature
                        'int32'  1 'a_cooler'; % Cooler
                        'int32'  1 'a_cooler_mode'; % Cooler mode
                        'int32'  1 'a_fan'; % Fan
                        ['uint',int2str(movie.data_depth)] [double(movie.width) double(movie.height)] 'IM'},...
                        'Repeat', N_frames_to_load);
                    
                    %keeping only the data
                    data=fmap.data;
                    clear fmap
                    
                    for i=1:N_frames_to_load
                        if data(i).magic ~= movie.magic,
                            ME = MException('moviereader:WrongMagicAtFrame',num2str(i));
                            throw(ME);
                        end
                        data(i).timestamp_sec = double(data(i).a_timestamp_sec) + double(data(i).a_timestamp_nsec)/10^9;
                        data(i).relative_timestamp_sec = data(i).timestamp_sec - movie.first_frame_timestamp_sec;
                    end
                    
                    %%%%% XIMEA %%%%%
                case CAMERA_TYPE_XIMEA
                    fmap=memmapfile([movie.Directory,movie.Filename],...
                        'Offset', offset,...
                        'Format', { 'uint32' 1 'magic';
                        'uint32' 1 'version'; %offset+4
                        'uint32' 1 'type'; %offset+8
                        'uint32' 1 'pixelmode'; %offset+12
                        'uint32' 1 'length_header'; %offset+16
                        'uint32' 1 'length_data'; %offset+20
                        'uint8' 100 'camera_name'; %offset+24 %this won't be readable, but I think 'chars' on MatLab are 2 bytes instead of 1 and I don't want to take any chance
                        'uint32' 1 'serial_number'; %offset+124
                        'uint64' 2 'timestamp'; %offset+128
                        'uint32' 1 'size_x_max'; %offset+144
                        'uint32' 1 'size_y_max'; %offset+148
                        'uint32' 1 'size_x'; %offset+152
                        'uint32' 1 'size_y'; %offset+156
                        'uint32' 1 'pos_x';
                        'uint32' 1 'pos_y';
                        'uint32' 1 'exposure';
                        'uint32' 1 'gain';  %this should be a 'float', but I think MatLab has an issue with them. They are 4 bytes as uint32, though.
                        'uint32' 1 'downsampling';
                        'uint32' 1 'downsampling_type';
                        'uint32' 1 'bad_pixel_correction';
                        'uint32' 1 'look_up_table';
                        'uint32' 1 'trigger';
                        'uint32' 1 'aeag';
                        'uint32' 1 'aeag_exposure_priority';    %should be a 'float'
                        'uint32' 1 'aeag_exposure_max_limit';
                        'uint32' 1 'aeag_gain_max_limit';      %should be 'float'
                        'uint32' 1 'aeag_average_intensity';
                        'uint32' 1 'hdr';
                        'uint32' 1 'hdr_t1';
                        'uint32' 1 'hdr_t2';
                        'uint32' 1 'hdr_t3';
                        'uint32' 1 'hdr_kneepoint1';
                        'uint32' 1 'hdr_kneepoint2';
                        
                        ['uint',int2str(movie.data_depth)] double([movie.width, movie.height]) 'IM'},...
                        'Repeat', N_frames_to_load);
                    
                    %keeping only the data
                    data=fmap.data;
                    clear fmap
                    
                    for i=1:N_frames_to_load
                        if data(i).magic ~= movie.magic,
                            ME = MException('moviereader:WrongMagicAtFrame',num2str(i));
                            throw(ME);
                        end
                        data(i).timestamp_sec = double(data(i).timestamp(1)) + double(data(i).timestamp(2))/10^9;
                        data(i).relative_timestamp_sec = data(i).timestamp_sec - movie.first_frame_timestamp_sec;
                    end
                    
            end
            
            
            % modifying the output so that it is back-compatible with
            % Davide's version
            
            
            timestamp = vertcat(data.relative_timestamp_sec);
            if N_frames_to_load > 3 %arbitrarily set
                movie.FrameRate = 1/mean(diff(timestamp));
            end
            
            
            if movie.data_depth == 8 || movie.data_depth == 16
                
                IM = [data.IM];
            elseif movie.data_depth == 12
                IM = vertcat(data.IM);
            end
            
            %             %check on magic word: if not TemI give an !!error!!
            %             wrong_magic = find(diff(vertcat(data.magic)),1,'first');
            %             if ~isempty(wrong_magic)
            %                 error(['File corrupted, wrong magic after ',num2str(wrong_magic+1)],' frames.');
            %             end
            
            
            % get field names, copy away the header. only makes sense if
            % all frames have same metadata
            headerfields = fieldnames(data(1));
            
            % copy the header fields from data structure
            for i=1:numel(headerfields)
                
                % skip the image data
                if any(strcmp(headerfields{i},{'IM','dummy'})), continue, end
                
                % copy the fields
                header.(headerfields{i}) = data(1).(headerfields{i});
                
            end %for
            
            clear data
            
            % reshaping the image vector
            IM=reshape(IM,movie.width,movie.height,N_frames_to_load);
            if size(IM,3)==1;
                IM=squeeze(IM);
            end
            
            IM = permute(IM, [2 1 3]);
            
            
            % e.g. to add output (remember you have to add them in
            % the output field of the function as well. Or you might want
            % to output a structure of arrrays (e.g. data_out with fields:
            % IM, timestamp, etc)
            % length_header = vertcat(data.length_header);
            % exposure = vertcat(data.exposure);
            % absolute_timestamp_sec = vertcat(data.timestamp_sec);
            
            %NB different cameras have different fields in "data". consider
            %moving the above lines inside the case stamements
            
        end
        
        
    end
    
    
end

