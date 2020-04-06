function [ Res, cccount ] = multiDDM_core( frame_stack, N_couple_frames_to_average,...
    Res)
%multiDDM_core Computational core for DDM calculation
%   Runs DDM on multiple tilings, as specified by Res. This function should
%   only really be called by the DDM_Analysis class.
%   Look at DDM_core in the supplementary information of the Original Work
%   listed below for DDM on full-frame.


%{
% Version 1.0
% © Luigi Feriani, Maya Juenet, Pietro Cicuta, 2018 (luigi.feriani@gmail.com) 
% 
% multiDDM_core.m is licensed under a Creative Commons 
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
%
% DDM algorithm based on 
% Cerbino, Trappe, 2008,  Physical review letters 100.18 (2008): 188102.
%
%}

%% input_check

if nargin < 2 || isempty(N_couple_frames_to_average)
    N_couple_frames_to_average = 100;
else
    if isa(N_couple_frames_to_average,'char') && ~strcmpi(N_couple_frames_to_average,'allframes')
        error('Input argument not recognised.');
    end
end

% check if gpu available
try
    gpuArray(1);
    flag_gpu = true;
catch
    flag_gpu=false;
end

%% general-purpose variables common to all boxes

[frame_size(1), frame_size(2), N_frames] = size(frame_stack);
max_tau = floor(N_frames/2);    %also equal to the number of lags

boxsize_vector = [Res.BoxSize];
N_boxsizes = numel(Res);


%% variables per box

% pixels per box
N_px_frame = boxsize_vector.^2;
fft2_norm_factor = 1./N_px_frame;

% save max_q once per boxsize
max_q = floor(boxsize_vector./2);
    
% number of boxes that fit in the field of view given BoxSize
N_Boxes_row = floor(frame_size(1) ./ boxsize_vector);
N_Boxes_col = floor(frame_size(2) ./ boxsize_vector);

row_span = N_Boxes_row .* boxsize_vector;
col_span = N_Boxes_col .* boxsize_vector;

row_offset = [Res.row_offset];
col_offset = [Res.col_offset];

% how many boxes per each boxsize?
N_boxes = N_Boxes_row .* N_Boxes_col;

% loop on boxsizes
for rc = N_boxsizes:-1:1

    % now save starting point of each box
    for j = N_Boxes_col(rc):-1:1
        for i = N_Boxes_row(rc):-1:1
            
            % index of box
            ii = sub2ind([N_Boxes_row(rc), N_Boxes_col(rc)], i, j);
            
            DRes(rc,1).Box(ii,1).row_start = Res(rc).row_offset + (i-1) * boxsize_vector(rc);
            DRes(rc).Box(ii).row_end   = Res(rc).row_offset + i * boxsize_vector(rc) - 1;
            
            DRes(rc).Box(ii).col_start = Res(rc).col_offset + (j-1) * boxsize_vector(rc);
            DRes(rc).Box(ii).col_end   = Res(rc).col_offset + j * boxsize_vector(rc) - 1;
            
            % initialise all the Iqtaus
            Res(rc).Box(ii).Iqtau = zeros(max_q(rc), max_tau, 'double');
            
            
        end %for i
    end %for j
    
    % now create and save the distance maps for fast radial average
    [DRes(rc).distance_map, DRes(rc).dist_counts] = distmapfun(boxsize_vector(rc));
end



%% actual DDM calculation
%(difference of frames, then |FFT|^2, average and radial average)


cccount = zeros(1,max_tau);

for tau = 1:max_tau
    
    % initialise the matrix for running averages
    for rc = 1:N_boxsizes
        DRes(rc).accum_abs_FT_diff_image = ...
            zeros([boxsize_vector(rc),boxsize_vector(rc),N_boxes(rc)], 'single');
    end
    
    
    ind_frames = choose_ind_frames(tau, N_couple_frames_to_average, N_frames);
    ccc = numel(ind_frames); % number of couples of frames at this tau
    
    for i = ind_frames % on-the-fly average of fft2 of difference images
        
        % take difference
        if flag_gpu
            tempdiff = gpuArray(single(frame_stack(:,:,i)) - single(frame_stack(:,:,i+tau)));
        else
            tempdiff = single(frame_stack(:,:,i)) - single(frame_stack(:,:,i+tau));
        end
        
        % now take the right fft of the right bit for each box
        for rc = 1:N_boxsizes
            
            tempdiff_chunk = tempdiff(row_offset(rc):row_offset(rc)+row_span(rc)-1,...
                col_offset(rc):col_offset(rc)+col_span(rc)-1);
            
            % reshape into 3d stack for fast fft2
            tempdiff_stack = reshape(permute(...
                reshape(tempdiff_chunk, row_span(rc), boxsize_vector(rc),[]), [2 1 3]),...
                boxsize_vector(rc), boxsize_vector(rc), []);
            
            % fft2 of difference image
            ft_tempdiff = fft2( tempdiff_stack ) * fft2_norm_factor(rc);
            
            
            % accumulate
            if flag_gpu
                DRes(rc).accum_abs_FT_diff_image = ...
                    DRes(rc).accum_abs_FT_diff_image +...
                    gather(real( ft_tempdiff ).^2 + imag( ft_tempdiff ).^2);
            else
                DRes(rc).accum_abs_FT_diff_image = ...
                    DRes(rc).accum_abs_FT_diff_image +...
                    real( ft_tempdiff ).^2 + imag( ft_tempdiff ).^2;
            end %if
            
        end %for rc
        
    end % for ind_frames
    
    
    % divide by number of couples of frames at this tau, take the
    % radial averages and write result in Iqtaus
    for rc = 1:N_boxsizes
        
        %average on initial times
        DRes(rc).averaged_abs_FT_diff_image = ...
            DRes(rc).accum_abs_FT_diff_image./ccc;
        
        for bc = 1:N_boxes(rc)
            
            % radial average (no need to store this
            oneD_power_spectrum = ...
                accumarray(DRes(rc).distance_map,...
                reshape(DRes(rc).averaged_abs_FT_diff_image(:,:,bc),[],1) )./ ...
                DRes(rc).dist_counts;	% radial average
            
            % fill each column of the output with the 1D power spectrum
            % just calculated
            Res(rc).Box(bc).Iqtau(:,tau) = ...
                oneD_power_spectrum(2:max_q(rc)+1);
            
        end
    end
    
    cccount(tau) = ccc;
    
    

end

end %multiDDM_core


function [distance_map, dist_counts] = distmapfun(frame_size)
%distance map for fast radial average

jj = repmat((1:frame_size),frame_size,1);
ii=jj';
cc = floor(frame_size/2)+1;
distance_map = fftshift(round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1); %if we fftshit here is only 1 time instead of Nframes/2
distance_map = distance_map(:);

dist_counts = accumarray(distance_map,ones(frame_size*frame_size,1));

end %function


function ind_frames = choose_ind_frames(tau, N_couple_frames_to_average, N_frames)
%choose_ind_frames returns the index of frames to use for subtracting

if isa(N_couple_frames_to_average,'char')
    
    ind_frames = 1:N_frames - tau;
    
else
    
    ind_frames_1 = randi(tau,1):tau:N_frames-tau;
    ind_frames_2 = randi(N_frames-tau,[1,N_couple_frames_to_average]);
    
    if numel(ind_frames_2)<numel(ind_frames_1)
        ind_frames = ind_frames_2;
    else
        ind_frames = ind_frames_1;
    end
    
    % maybe change the last lines to:
%     ind_frames = randi(tau,1):tau:N_frames-tau;
%     n = numel(ind_frames);
%     if numel(ind_frames) > N_couple_frames_to_average % if too many
%         %pick N_couple_frames_to_average at random
%         ind_frames = ind_frames(randi(n,[1,N_couple_frames_to_average]));
%     end

end

end %function
