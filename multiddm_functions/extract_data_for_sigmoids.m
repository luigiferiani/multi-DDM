function [AccumData] = extract_data_for_sigmoids(filelist, boxsizes_vector, q_limits_1oum, flag_recalculate_goodboxes)
%extract_data_for_sigmoids loops the files in the structure (that comes
%from dir) filelist extracting the relevant data (decay and frequency) for
%all the boxsizes in boxsizes_vector for the modes within q_limits_1omum

% range of q. This depends on objective(s) used in a sigmoid
% q_limits_1oum = [1.95 2.75]; % have to use this for comparing, 20X, 40X and 60X

%{
% Version 1.0
% © Luigi Feriani 2019 (luigi.feriani@gmail.com) 
% 
% extract_data_for_sigmoids.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% 
% Original work:
% 
% Chioccioli, M.*, Feriani, L.*, Kotar, J., Bratcher, P. E.**, Cicuta, P.**, Nature Communications 2019
% "Phenotyping ciliary dynamics and coordination in response to CFTR-modulators 
% in Cystic Fibrosis respiratory epithelial cells"
%}

%% input check

% handle the case in which the filelist is empty. in this case
% boxsizes_vector cannot be empty!
if isempty(filelist)
    
    N_BoxSizes = numel(boxsizes_vector);
    
    % cell arrays (N_BoxSizes x 1)
    AccumData.Frequency_Hz      = cell(N_BoxSizes,1);
    AccumData.Damping_Hz        = cell(N_BoxSizes,1);
    AccumData.Damping_1ocycles  = cell(N_BoxSizes,1);
    
    % matrices (N_BoxSizes x number of files)
    AccumData.good_fits_fraction = [];
    
    % array (N_BoxSizes x 1)
    AccumData.BoxSizes_px = boxsizes_vector(:);
    AccumData.px2mum = [];
    AccumData.BoxSizes_um = [];
    
%     warning('Empty file list, empty structure returned.');
    return
end 

if nargin < 4 || isempty(flag_recalculate_goodboxes)
    flag_recalculate_goodboxes = false;
end

if nargin < 3 || isempty(q_limits_1oum)
    warning('q limits empty - using default one. This works for Grasshopper, 40x, 1920x1024.');
    q_limits_1oum = [1.95 2.75];
end

if nargin < 2 || isempty(boxsizes_vector)
    warning('Box size vectors not specified, reading from the 1st file in the list');
    temp = load(fullfile(filelist(1).folder,filelist(1).name));
    boxsizes_vector =  sort([temp.cilia.Results(:).BoxSize]);
%     boxsizes_vector =  [temp.cilia.Results(:).BoxSize];
end



%% Initialise accumulant structures

N_BoxSizes = numel(boxsizes_vector);

% cell arrays, one growing vector per each box size

% create variables to put together all results from different files
accum_frequency_Hz = cell(N_BoxSizes,1);
accum_damping_Hz = cell(N_BoxSizes,1);
accum_damping_1ocycles = cell(N_BoxSizes,1);
accum_file_ind = cell(N_BoxSizes,1);
accum_box_ind = cell(N_BoxSizes,1);

% create variables that tells us the fraction of good fits per video per
% boxsize, it's a number of boxsizes-by-number of files matrix
good_fits_fraction = nan(N_BoxSizes,numel(filelist));
good_boxes_fraction = nan(N_BoxSizes,numel(filelist));
movement_boxes_fraction = nan(N_BoxSizes,numel(filelist));

%% -------------- PUT DATA TOGETHER ---------------- %%

% display file list
try
    disp(vertcat(filelist.name));
catch
end

% this is for control, to check that we aren't putting together videos with
% different px2mum as we rely on size in pixels for categorising the
% plots
px2mum_vector = zeros(size(filelist));

% for loop on files
for fc = 1:numel(filelist)
    
    % output string
    ostr_counter = sprintf('%.2d/%.2d',fc,numel(filelist));
    fprintf(ostr_counter);
    
    load(fullfile(filelist(fc).folder,filelist(fc).name));
    
    % make this happen with a flag!
    if flag_recalculate_goodboxes
        cilia.gather_results(true);
    end %if
    
    % save px2mum
    px2mum_vector(fc) = cilia.px2mum;
    if fc >= 2 && px2mum_vector(fc) ~= px2mum_vector(fc-1)
        error('Not all the videos had the same magnification. Please restrict the file list and retry.');
    end
    
    
    % loop on the boxsizes_vector
    
    for bsc = 1:numel(boxsizes_vector) %bsc = boxsize_counter
        
        % current box size
        bs = boxsizes_vector(bsc);
        
        % index of current box size in the results structure
        bi = find( [cilia.Results.BoxSize] == bs, 1);
        
        % can be empty, in case move on
        if isempty(bi)
            continue
        end %if
        
        % find the mode range corresponding to the q range
        q_limits_1opx = q_limits_1oum * cilia.px2mum;
        idx_q_inrange = isinrange(cilia.Results(bi).qVec, q_limits_1opx );
        
        % exclude bad boxes (threshold of std_fs) first
        % these matrices are (number of modes x number of good boxes) big
        temp_Frequency_Hz_mat = horzcat(cilia.Results(bi).Box(cilia.Results(bi).ind_good_boxes).Frequency).*cilia.FrameRate; % in Hz
        temp_Damping_Hz_mat = horzcat(cilia.Results(bi).Box(cilia.Results(bi).ind_good_boxes).Damping).*cilia.FrameRate; %in Hz
        temp_Amplitude_mat = horzcat(cilia.Results(bi).Box(cilia.Results(bi).ind_good_boxes).Amplitude);
        temp_box_ind_mat = repmat(1:nnz(cilia.Results(bi).ind_good_boxes), size(temp_Frequency_Hz_mat,1),1);
        
        % this calculates all the sum of squares from fit for the good
        % boxes
        temp_norm_sumsquares_mat = zeros(size(temp_Frequency_Hz_mat));
        counter = 0;
        for bc = find(cilia.Results(bi).ind_good_boxes(:))'
            counter = counter + 1;
            max_q_fitted = cilia.Results(bi).BoxSize/4;
            % using "norm" data should ensure we get rid of possible
            % mistakes in normalising Iqtau at different boxsizes
            norm_data = (cilia.Results(bi).Box(bc).Iqtau(1:max_q_fitted,1:cilia.max_tau_fitted) - repmat( cilia.Results(bi).Box(bc).Offset,1,cilia.max_tau_fitted) ) ./ repmat( cilia.Results(bi).Box(bc).Amplitude,1,cilia.max_tau_fitted);
            norm_fit = (1 - exp(cos( 2*pi* cilia.Results(bi).Box(bc).Frequency * (1:cilia.max_tau_fitted)  ))) .* exp( - cilia.Results(bi).Box(bc).Damping * (1:cilia.max_tau_fitted));
            norm_sumsquares = sum(abs(norm_data - norm_fit),2);
            temp_norm_sumsquares_mat(:,counter) = norm_sumsquares;
        end %for
        
        
        % restrict the temp matrices to the selected q range => these are
        % now number of modes in range x number of good boxes
        try
        temp_Frequency_Hz_mat = temp_Frequency_Hz_mat(idx_q_inrange,:);
        temp_Damping_Hz_mat = temp_Damping_Hz_mat(idx_q_inrange,:);
        temp_Amplitude_mat = temp_Amplitude_mat(idx_q_inrange,:);
        temp_box_ind_mat = temp_box_ind_mat(idx_q_inrange,:);
        temp_norm_sumsquares_mat = temp_norm_sumsquares_mat(idx_q_inrange,:);
        catch
            filelist(fc).name
            idx_q_inrange
            temp_Frequency_Hz_mat
            nnz(cilia.Results(bi).ind_good_boxes(:))
            fprintf('\n');
        end
        
        % control to throw stuff away
        bad_fits = false(size(temp_Frequency_Hz_mat));
        
        % controls on frequency
        bad_fits(temp_Frequency_Hz_mat <=   1)  = true;    %exclude where frequency too low (aka sth went wrong in the fit or the sample is still)
        bad_fits(temp_Frequency_Hz_mat >=  30)  = true;    %exclude where frequency too high (aka sth went wrong in the fit)
        bad_fits(isnan(temp_Frequency_Hz_mat))  = true;    %exclude where frequency is nan (aka sth went wrong in the fit)
        
        % controls on damping
%         bad_fits(temp_Damping_mat   > 2*temp_Frequency_mat)  = true;    %exclude where decay time is less than a semiperiod (aka damping is more than 2*frequency) because in this case the frequency is probably wrong
        bad_fits(temp_Damping_Hz_mat   > 3*temp_Frequency_Hz_mat)  = true;    %exclude where decay time is less than a 3rd of a period (aka damping is more than 3*frequency) because in this case the frequency is probably wrong
        bad_fits(temp_Damping_Hz_mat   == 0)  = true;    %exclude where decay time is 0 (aka sth went wrong in the fit)
        bad_fits(temp_Damping_Hz_mat < 1/(10*cilia.NumberOfFrames/cilia.FrameRate)) = true;  %exclude where decay time is more than 10*length of video
        bad_fits(isnan(temp_Damping_Hz_mat))  = true;            %exclude where damping is nan (aka sth went wrong in the fit)
        
        % controls on amplitude
        bad_fits(temp_Amplitude_mat <= eps)  = true;    %exclude where signal too low
        bad_fits(isnan(temp_Amplitude_mat))  = true;    %exclude where amplitude is nan (aka sth went wrong in the fit)
        
        
        
        % controls on sum of squares
        %         % this worked on cedar's data
        %         bad_fits(temp_norm_sumsquares_mat > 30) = true; %excludes fit with too high sum of distances data-fit
        % sum of squares are roughly lognormal. Find center, sigma, and exclude fits that are out of 2 sigma away on the right hand side
        mm = mean(log10(temp_norm_sumsquares_mat(~bad_fits)));
        ss = std(log10(temp_norm_sumsquares_mat(~bad_fits)));
        bad_fits(log10(temp_norm_sumsquares_mat) > mm+2*ss) = true;
        
        % change notation
        good_fits = ~bad_fits;
        
        % calculates how many good fits in the q range selected for this
        % file and boxsize
        good_fits_fraction(bsc,fc) = sum(good_fits(:)) / numel(good_fits);

        % how many good boxes? fraction of good boxes out of the initial
        % good boxes times the initial fraction of good boxes. So basically
        % the good boxes here (any(good_fits)) divided by the total number
        % of boxes
        good_boxes_fraction(bsc,fc) = sum(any(good_fits,1)) / numel(cilia.Results(bi).ind_good_boxes);
        movement_boxes_fraction(bsc,fc) = nnz(cilia.Results(bi).ind_good_boxes) / numel(cilia.Results(bi).ind_good_boxes);
        
        % grow the accumulant vectors:
        % put the column vector with good data together with the other good
        % data (good fit, good range) of same boxsize coming from the other videos
        accum_frequency_Hz{bsc} = vertcat(accum_frequency_Hz{bsc},...
            reshape(temp_Frequency_Hz_mat(good_fits(:)),[],1) );        
        accum_damping_Hz{bsc}   = vertcat(accum_damping_Hz{bsc},...
            reshape(temp_Damping_Hz_mat(good_fits(:)),[],1)   );
        accum_damping_1ocycles{bsc} = vertcat(accum_damping_1ocycles{bsc}, ...
            reshape(temp_Damping_Hz_mat(good_fits(:))./temp_Frequency_Hz_mat(good_fits(:)),[],1) );
        accum_file_ind{bsc} = vertcat(accum_file_ind{bsc},...
            fc .* ones(nnz(good_fits),1));
        accum_box_ind{bsc} = vertcat(accum_box_ind{bsc}, reshape(temp_box_ind_mat(good_fits),[],1));
        
    end %for

    % clean up
    clear cilia
    fprintf(repmat('\b',size(ostr_counter)));
    
end %for

%% Create structure for output 


% cell arrays (N_BoxSizes x 1)
AccumData.Frequency_Hz_fornorm      = accum_frequency_Hz;  %these are the ones used for normalising the damping
AccumData.Damping_Hz        = accum_damping_Hz;
AccumData.Damping_1ocycles  = accum_damping_1ocycles;
AccumData.file_ind          = accum_file_ind;        
AccumData.box_ind           = accum_box_ind; 

% the frequency for histograms instead should be averaged within each box
for bsc = 1:numel(boxsizes_vector)
    try
    temp_boxav_freq = accumarray({accum_file_ind{bsc}, accum_box_ind{bsc}}, accum_frequency_Hz{bsc},[],[], NaN) ./ ...
        accumarray({accum_file_ind{bsc}, accum_box_ind{bsc}}, ones(size(accum_frequency_Hz{bsc})),[],[], NaN);
    catch
        % put debugger here
    end
    temp_boxav_freq = temp_boxav_freq(~isnan(temp_boxav_freq));
    AccumData.Frequency_Hz{bsc,1} = temp_boxav_freq(:);
end %for

% matrices (N_BoxSizes x number of files)
AccumData.good_fits_fraction = good_fits_fraction;
AccumData.good_boxes_fraction = good_boxes_fraction;
AccumData.movement_boxes_fraction = movement_boxes_fraction;

% array (N_BoxSizes x 1)
AccumData.BoxSizes_px = boxsizes_vector(:);
AccumData.px2mum = unique(px2mum_vector);
AccumData.BoxSizes_um = AccumData.BoxSizes_px .* AccumData.px2mum;

% q limits
AccumData.q_limits_1oum = q_limits_1oum;
AccumData.q_limits_1opx = q_limits_1opx;

end %function

