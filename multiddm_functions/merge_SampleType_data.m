function [MergedData] = merge_SampleType_data(SampleType, ind_ST, ind_TP, donors_str, inserts_str, positions_str)
%merge_SampleType_data Takes the structure SampleType and indices to
%which set of data to merge, and produces a merged structure for plots
%   ind_ST          is an integer <= numel(SampleType)
%   ind_TP          is an integer <= numel(SampleType(i).TimePoints)
%   donors_str      is a char array or a cell array of char arrays
%   inserts_str     is a char array or a cell array of char arrays
%   positions_str   is a char array or a cell array of char arrays

%{
% Version 1.0
% © Luigi Feriani 2019 (luigi.feriani@gmail.com) 
% 
% merge_SampleType_data.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% 
% Original work:
% 
% Chioccioli, M.*, Feriani, L.*, Kotar, J., Bratcher, P. E.**, Cicuta, P.**, Nature Communications 2019
% "Phenotyping ciliary dynamics and coordination in response to CFTR-modulators 
% in Cystic Fibrosis respiratory epithelial cells"
%}

%% input check

if numel(ind_ST) > 1
    MergedData = [];
    error('Only one sample type!');
end

if numel(ind_TP) > 1
    MergedData = [];
    error('Only one timepoint!');
end

if ~sum( ~arrayfun(@(i)isempty(SampleType(ind_ST).TimePoint(ind_TP).Data(i).FileList),...
        1:numel(SampleType(ind_ST).TimePoint(ind_TP).Data) ) ) %this is true if they are all empty (~isempty is 0, its sum 0)
    disp('No data for this sample type at this timepoint.');
    MergedData = [];
    return
end

% force strings to become cells
donors_str  = cellstr(donors_str);
inserts_str = cellstr(inserts_str);
positions_str = cellstr(positions_str);

flag_debugging = false;

%% make a copy of useless info (for my convenience while coding)

% sample type
sampletype_str = SampleType(ind_ST).Str;

% timepoint
timepoint_str = SampleType(ind_ST).TimePoint(ind_TP).Str;

% data structure
TMP = SampleType(ind_ST).TimePoint(ind_TP).Data;    %this has N_averageable_categories entries


%% select the appropriate data to merge as per input

% which entries of Data do we merge?

% entries of data full in the first place
idx_present_data = ~arrayfun(@(i)isempty(TMP(i).FileList),(1:numel(TMP))');

% compare donors
idx_merge_donors = false(size(TMP));
for i = 1:numel(donors_str)
    idx_merge_donors = idx_merge_donors | strcmp({TMP.Donor}',donors_str{i});
end

% compare inserts
idx_merge_inserts = false(size(TMP));
for i = 1:numel(inserts_str)
    idx_merge_inserts = idx_merge_inserts | strcmp({TMP.Insert}',inserts_str{i});
end

% compare positions
idx_merge_positions = false(size(TMP));
for i = 1:numel(positions_str)
    idx_merge_positions = idx_merge_positions | strcmp({TMP.Position}',positions_str{i});
end

% now combine
idx_merge = idx_present_data & idx_merge_donors & idx_merge_inserts & idx_merge_positions;
if sum(idx_merge) == 0
    MergedData = []
    return
end

%% merge the appropriate quantities

% maybe I can actually vectorise this, but not sure how to

% the order of box sizes should be the same! but i can put a check just in
% case

% this select just the averageable entries of TMP
AccumDataMerge = vertcat(TMP(idx_merge).AccumData);
N_merged = sum(idx_merge);

% check on boxsizes
if any(any(diff(horzcat(AccumDataMerge.BoxSizes_px),1,2)))
    error('The box sizes are not in the same order!');
else %then just take on of them
    % find a 
    window_area_um2 = AccumDataMerge(1).BoxSizes_um.^2; % will need this for x axis
end %if

% merge quantities
for bsc = 1:numel(AccumDataMerge(1).BoxSizes_px)
    
    Frequency_Hz{bsc,1} = vertcat(cell2mat(arrayfun(@(i)AccumDataMerge(i).Frequency_Hz{bsc},(1:N_merged)',...
        'UniformOutput',false)));
    
    Damping_Hz{bsc,1} = vertcat(cell2mat(arrayfun(@(i)AccumDataMerge(i).Damping_Hz{bsc},(1:N_merged)',...
        'UniformOutput',false)));
    
    Damping_1ocycles{bsc,1} = vertcat(cell2mat(arrayfun(@(i)AccumDataMerge(i).Damping_1ocycles{bsc},(1:N_merged)',...
        'UniformOutput',false)));
    
end %for bsc
    
% this is N_Boxsizes by number of files the data come from, i.e. numel(vertcat(TMP(idx_merge).FileList))
good_fits_fraction      = horzcat(AccumDataMerge.good_fits_fraction);
good_boxes_fraction     = horzcat(AccumDataMerge.good_boxes_fraction);
movement_boxes_fraction = horzcat(AccumDataMerge.movement_boxes_fraction);



if flag_debugging

horzcat(AccumDataMerge.Frequency_Hz)   
Frequency_Hz

horzcat(AccumDataMerge.Damping_Hz)   
Damping_Hz

horzcat(AccumDataMerge.Damping_1ocycles)  
Damping_1ocycles

numel(vertcat(TMP(idx_merge).FileList))
size(good_fits_fraction)

end %if


%% calculate the median damping for clean sigmoid plottinf

% median of not-normalised sigmoid
med_Damping_Hz = cellfun(@median,Damping_Hz);
ler_Damping_Hz = + med_Damping_Hz - cellfun(@(x)prctile(x,25), Damping_Hz);
uer_Damping_Hz = - med_Damping_Hz + cellfun(@(x)prctile(x,75), Damping_Hz);

% median of normalised sigmoid
med_Damping_1ocycles = cellfun(@median,Damping_1ocycles);
ler_Damping_1ocycles = + med_Damping_1ocycles - cellfun(@(x)prctile(x,25), Damping_1ocycles);
uer_Damping_1ocycles = - med_Damping_1ocycles + cellfun(@(x)prctile(x,75), Damping_1ocycles);


%% fit the median points

fo = fitoptions('Method','NonLinearLeastSquares','StartPoint', [ 9 1e3 0]);
ft = 'a/(1+b./x)+c';

fo2 = fitoptions('Method','NonLinearLeastSquares',...
    'StartPoint', [ max(med_Damping_Hz) min(med_Damping_Hz) 3.5],...
    'Lower',[0 0 0]);
ft2 = 'a/(1+10.^(mu-x))+c';

% med_Damping_* can contain nans, if any of the entries of Damping_Hz was
% empty. So have to account for that.


try
    inn = ~isnan(med_Damping_Hz); % idx not nans
    [Damping_Hz_fit_out,       Damping_Hz_gof]       = fit(window_area_um2(inn), med_Damping_Hz(inn),        ft,  fo);
    [Damping_Hz_fit_out2,      Damping_Hz_gof2]      = fit(log10(window_area_um2(inn)), med_Damping_Hz(inn), ft2, fo2);
    [Damping_1ocycles_fit_out, Damping_1ocycles_gof] = fit(window_area_um2(inn), med_Damping_1ocycles(inn),  ft,  fo);
catch
    sampletype_str
    timepoint_str
    unique({TMP(idx_merge).Donor})
    unique({TMP(idx_merge).Insert})
    unique({TMP(idx_merge).Position})
    window_area_um2
    med_Damping_Hz
    med_Damping_1ocycles
    return
end


%% prepare data for output

% strings
MergedData.sampletype_str       = sampletype_str;
MergedData.timepoint_str        = timepoint_str;
MergedData.donors_str           = unique({TMP(idx_merge).Donor});       % only the ones actually merged
MergedData.inserts_str          = unique({TMP(idx_merge).Insert});      % only the ones actually merged
MergedData.positions_str        = unique({TMP(idx_merge).Position});   % only the ones actually merged

% freq
MergedData.Frequency_Hz         = Frequency_Hz;

% sigmoids
MergedData.window_area_um2          = window_area_um2;

MergedData.Damping_Hz               = Damping_Hz;
MergedData.med_Damping_Hz           = med_Damping_Hz;
MergedData.ler_Damping_Hz           = ler_Damping_Hz;
MergedData.uer_Damping_Hz           = uer_Damping_Hz;
MergedData.Damping_Hz_fit_out       = Damping_Hz_fit_out;
MergedData.Damping_Hz_gof           = Damping_Hz_gof;
MergedData.Damping_Hz_fit_out2      = Damping_Hz_fit_out2;
MergedData.Damping_Hz_gof2          = Damping_Hz_gof2;

MergedData.Damping_1ocycles         = Damping_1ocycles;
MergedData.med_Damping_1ocycles     = med_Damping_1ocycles;
MergedData.ler_Damping_1ocycles     = ler_Damping_1ocycles;
MergedData.uer_Damping_1ocycles     = uer_Damping_1ocycles;
MergedData.Damping_1ocycles_fit_out = Damping_1ocycles_fit_out;
MergedData.Damping_1ocycles_gof     = Damping_1ocycles_gof;

% checks
MergedData.good_fits_fraction       = good_fits_fraction;
MergedData.good_boxes_fraction      = good_boxes_fraction;
MergedData.movement_boxes_fraction  = movement_boxes_fraction;



end %function

















