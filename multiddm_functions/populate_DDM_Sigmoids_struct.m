function [SampleType] = populate_DDM_Sigmoids_struct( analysis_folder, ...
    imaging_string, sampletypes, timepoints, donors, inserts, positions,  ...
    boxsizes_vector, q_limits_1oum, flag_dryrun, flag_recalculate_goodboxes)
%populate_DDM_Sigmoids_struct creates a structure array based on SampleTypes,
% TimePoints, loops through the files in AnalysisFolder, and categorises
% them into the structure array using the proper Donors, Inserts, Positions
% labels. This function will not have a GUI, but will probably come after
% the parameters have been set via a GUI
%
%   AnalysisFolder is a cell array of chars, as data can be in more than
%   one folder
%
%   ImagingString is a string that has to be shared by *all* files. Usually
%   contains magnification and type of imaging (e.g. BF). Can just be '40X', or '40X_BF' and it's
%   still ok
%
%   SampleTypes can be a cell array of chars, (e.g. {'CF', 'CFm', 'N',
%   'Nm'} etc
%
%   TimePoints can be a cell array of chars. If there is no use for it,
%   then just put the day or leave it empty
%
%   Donors, Inserts, Positions can be cell arrays, no need to fill them
%   all, can be left empty
%

%{
% Version 1.0
% © Luigi Feriani 2019 (luigi.feriani@gmail.com) 
% 
% populate_DDM_Sigmoids_struct.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% 
% Original work:
% 
% Chioccioli, M.*, Feriani, L.*, Kotar, J., Bratcher, P. E.**, Cicuta, P.**, Nature Communications 2019
% "Phenotyping ciliary dynamics and coordination in response to CFTR-modulators 
% in Cystic Fibrosis respiratory epithelial cells"
%}


%% input check

if nargin < 11 || isempty(flag_recalculate_goodboxes)
    flag_recalculate_goodboxes = false;
end

if nargin < 10 || isempty(flag_dryrun)
    flag_dryrun = false;
end

if nargin < 9 || isempty(q_limits_1oum)
    warning('q limits empty - using default one. This works for Grasshopper, 20x, 40x or 60x, 1920x1024.');
    q_limits_1oum = [1.95 2.75];
end

if nargin < 8 || isempty(boxsizes_vector)
%     warning('Box size vectors not specified, using the default "quick" one');
%     boxsizes_vector =  [32 64 128 256 512 1024];
    boxsizes_vector = [];
end

if nargin < 7 || isempty(positions)
    positions = '';
end %if
if iscellstr(positions) || ischar(positions)
    positions = cellstr(positions);
    positions = positions(:);
    Npos = numel(positions);
else
    error('Positions has to be a string or a cell array of strings.');
end %if


if nargin < 6 || isempty(inserts)
    inserts = '';
end %if
if iscellstr(inserts) || ischar(inserts)
    inserts = cellstr(inserts);
    inserts = inserts(:);
    Nins = numel(inserts);
else
    error('Inserts has to be a string or a cell array of strings.');
end %if

if nargin < 5 || isempty(donors)
    donors = '';
end %if
if iscellstr(donors) || ischar(donors)
    donors = cellstr(donors);
    donors = donors(:);
    Ndon = numel(donors);
else
    error('Donors has to be a string or a cell array of strings.');
end %if

if nargin < 4 || isempty(timepoints)
    timepoints = '';
end %if
if iscellstr(timepoints) || ischar(timepoints)
    timepoints = cellstr(timepoints);
    timepoints = timepoints(:);
    Ntpt = numel(timepoints);
else
    error('TimePoints has to be a string or a cell array of strings.');
end %if

if nargin < 3 || isempty(sampletypes)
    error('You need to provide SampleTypes')
else
    if iscellstr(sampletypes) || ischar(sampletypes)
        sampletypes = cellstr(sampletypes);
        sampletypes = sampletypes(:);
        Nstp = numel(sampletypes);
    else
        error('SampleTypes has to be a string or a cell array of strings.');
    end %if
    
end %if

if nargin < 2 || isempty(imaging_string)
    imaging_string = '';
else
    if ~ischar(imaging_string)
        error('ImagingString has to be a string.');
    end %if
end %if

if nargin < 1 || isempty(analysis_folder)
    analysis_folder = pwd;
    Nafd = 1;
end %if
if iscellstr(analysis_folder) || ischar(analysis_folder)
    analysis_folder = cellstr(analysis_folder);
    analysis_folder = analysis_folder(:);
    Nafd = numel(analysis_folder);
else
    error('AnalysisFolders has to be a string or a cell array of strings.');
end %if


flag_debugging = true;

%% prepare categories

% how many averageable categories? if I can average on donors, inserts,
% positions then I have donors * inserts * positions
Navcat = Ndon * Nins * Npos;

% create a Navcat * 3 cell array with the combination
cc = 1;
for dnc = 1:Ndon
    for inc = 1:Nins
        for  psc = 1:Npos
            AvCat{cc,1} = donors{dnc};
            AvCat{cc,2} = inserts{inc};
            AvCat{cc,3} = positions{psc};
            cc = cc+1;
        end %for psc
    end %for inc
end %for dnc


%% create a full list of files

for afc = 1:Nafd
    
    temp(afc).fl = dir(fullfile(analysis_folder{afc},[imaging_string,'*.mat']));
    
end %for
    
% unsorted list of all files
init_fl = vertcat(temp.fl); % initial_filelist
    

%% associate the file list with lists of categories

init_fl_ST = match_filelist_category(init_fl, sampletypes);
init_fl_TP = match_filelist_category(init_fl, timepoints);
init_fl_D = match_filelist_category(init_fl, donors);
init_fl_I = match_filelist_category(init_fl, inserts);
init_fl_P = match_filelist_category(init_fl, positions);

% ouptut for debugging
if flag_debugging
for i = 1:numel(init_fl)
    fprintf('\n%s\t\t%s\t%s\t%s\t%s\t%s',init_fl(i).name, init_fl_ST{i}, init_fl_TP{i}, init_fl_D{i}, init_fl_I{i}, init_fl_P{i})
end
end %if
fprintf('\n');

%% create the structure

for stc = Nstp:-1:1

    SampleType(stc,1).Str = sampletypes{stc};
    
    % filenames that match the current sampletype
%     idx_fl_ST = strcmp(init_fl_ST, SampleType(stc,1).Str);
    idx_fl_ST = strcmp_wildcard(init_fl_ST, SampleType(stc,1).Str);

    for tpc = Ntpt:-1:1

        SampleType(stc).TimePoint(tpc,1).Str = timepoints{tpc};
    
        % filenames that match the current timepoints
%         idx_fl_TP = strcmp_wildcard(init_fl_TP,
%         SampleType(stc).TimePoint(tpc,1).Str); %this doesn't work with
%         day10, afterwash_day10
        idx_fl_TP = strcmp(init_fl_TP, SampleType(stc).TimePoint(tpc,1).Str);

            
        % now I hoop on the combinations of averageable categories
        for cc = 1:Navcat
            
            % reset
            clearvars Data
            
            % Averageable categories
            Data.Donor      = AvCat{cc,1};
            Data.Insert     = AvCat{cc,2};
            Data.Position   = AvCat{cc,3};
            
            % filenames matching the averageable categories
            idx_fl_D = strcmp(init_fl_D, Data.Donor);
            idx_fl_I = strcmp(init_fl_I, Data.Insert);
            idx_fl_P = strcmp(init_fl_P, Data.Position);

            
            % finally write only the right files            
            Data.FileList  = init_fl(idx_fl_ST & idx_fl_TP & idx_fl_D & idx_fl_I & idx_fl_P);
            
            % actually go and find data. Takes time.
            if ~flag_dryrun
                Data.AccumData =  extract_data_for_sigmoids(Data.FileList, boxsizes_vector,...
                    q_limits_1oum, flag_recalculate_goodboxes);
            end %if
            
            % now store in output structure
            SampleType(stc).TimePoint(tpc).Data(cc,1) = Data;
            
        end %for cc
    end %for tpc
end %for stc

try
    structstruct(SampleType);
catch
end









end %function




function [filelist_matches] = match_filelist_category(filelist, categories)
%match_filelist_category creates a cell array of size(filelist) with the
%match of regexp(filelist(i).name, categories,'match','once'), after
%translating categories into regexps

% input check
if ischar(categories)
    categories = cellstr(categories);
end
categories = categories(:);

% translate by replacing the standard * with the regexp .*
categories = regexptranslate('wildcard',categories);

% initialise
filelist_matches = cell(size(filelist)); %initial_filelist_SampleType


% for loop on files
for fc = 1:numel(filelist)
    
    temp_match = regexp(filelist(fc).name,categories,'match','once');
    
    % three cases: either there's no matches, or there's one, or there's
    % more than one. temp_out will always have Nstp elements, they might be
    % empty though
    
%     N_matches = sum(arrayfun(@(i)~isempty(temp_match{i}),1:Nstp));
    N_matches = sum(~cellfun(@isempty,temp_match));
    
    if N_matches == 0 % no matches found. This shouldn't happen
        
        filelist_matches{fc} = '';
        
    elseif N_matches == 1 % only one match found
        
        filelist_matches{fc} = temp_match{~cellfun(@isempty,temp_match)};
        
    elseif N_matches >=2 
        % more than a match found. Unless the input contained two disjoint
        % strings that were both found in the filename, then the longer
        % match is the good one
        
        % measure length of matches
        temp_match_length = cellfun(@length,temp_match);       
        [~,ind_match] = max(temp_match_length);
        filelist_matches{fc} = temp_match{ind_match};
        
    end %if
    
end %for fc

end %function




function [tf] = strcmp_wildcard(input, pattern)
%strcmp_wildcard looks for matches of input with pattern. pattern can
%contain wildcards. Returns true when pattern is found in input, false
%otherwise. input can be a cell array

% translate by replacing the standard * with the regexp .*
pattern = regexptranslate('wildcard',pattern);

% search
matches = regexp(input, pattern, 'match','once');

% find successful matches
tf = ~cellfun(@isempty, matches);


end %function














