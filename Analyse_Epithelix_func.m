function [] = Analyse_Epithelix_func(data_folder, analysis_folder, files_string, boxes_mode)
%Analyse_Epithelix_func allows the user to select a data folder, analysis
%folder, files identification string, and sizes of the DDM windows, and
%runs instances of the DDM_Analysis on all files thus selected.
%
%   Analyse_Epithelix_func() opens a GUI asking the user to input the
%   variables needed
%
%   Analyse_Epithelix_func(data_folder, analysis_folder, files_string,
%   boxes_mode) starts the analysis without the GUI unless:
%   a. one of the inputs is empty,
%   b. files_string = '*.movie', or
%   c. boxes_mode   = 'custom'.
%
%   boxes_mode has to either be 'quick', 'sigmoids', or 'custom'.


%{
% Version 1.1
% © Luigi Feriani 2019 (luigi.feriani@gmail.com) 
% 
% Analyse_Epithelix_func.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% 
% Original work
% Feriani, L., et al., Biohpysical Journal 2017
% "Assessing the Collective Dynamics of Motile Cilia in Cultures of Human
% Airway Cells"
%
%}

%% input check and settings

splash;

% data folder
if nargin < 1 || isempty(data_folder)
    data_folder = [];
    fprintf('\nData folder not passed as input, manual selection.');
else
    % check if local path
    if ~contains(data_folder,'/') && ~contains(data_folder,'\') % then relative
        data_folder = fullfile(pwd,data_folder);
    end
end %if

% analysis folder
if nargin < 2 || isempty(analysis_folder)
    analysis_folder = [];
else
    % check if local path
    if ~contains(analysis_folder,'/') && ~contains(analysis_folder,'\') % then relative
        analysis_folder = fullfile(pwd,analysis_folder);
    end
end %if

% files string
if nargin < 3 || isempty(files_string)
    files_string = '*.movie';
end %if

% boxsizes vector
if nargin < 4 || isempty(boxes_mode) || ~ischar(boxes_mode)
    boxes_mode = 'custom';
end
boxsizes_vector = get_boxsizes_vector(boxes_mode);

if isempty(data_folder) || isempty(analysis_folder) || strcmp(boxes_mode,'custom') || strcmp(files_string,'*.movie')
    [data_folder, analysis_folder, files_string, boxsizes_vector ] =...
        GUI_settings(data_folder, analysis_folder, files_string, boxes_mode);
end %if


%% summary output

Summary = [  sprintf('Summary:\n'),...
    sprintf( 'Computer        = %s\n',getenv('computername')),...
    sprintf( 'Data folder     = %s\n',data_folder),...
    sprintf( 'Analysis folder = %s\n',analysis_folder),...
    [sprintf('Box sizes       = '),  sprintf('%d  ',boxsizes_vector)], sprintf('\n'),...
    sprintf( 'Files string    = %s\n',files_string) ];

fprintf('%s',Summary);


%% prepare for execution

if ~isdir(analysis_folder)
    mkdir(analysis_folder);
end

filelist = dir(fullfile(data_folder,files_string));


%% execution (do not modify!!)

tic;
timeelapsed = zeros(numel(filelist),1);
startfrom = 1;
stopat = numel(filelist);

try
    for i = startfrom:stopat
        
        % create fullfilename
        filename = fullfile(data_folder, filelist(i).name);
        [~,savename,~] = fileparts(filelist(i).name);
        savename = fullfile(analysis_folder,[savename,'.mat']);
        
        cprintf('*[0 0 .4]','\n%d/%d',i,stopat)
        cprintf('*[0 0 .4]','\n%s',filename)
        cprintf('*[0 0 .4]','\n%s\n',savename)
        
        if isempty(dir(savename))
            cilia = DDM_Analysis(filename);
            cilia.set_temperature(1);
            cilia.N_couple_frames_to_average = 200;
            cilia.VariableBoxSize_Analysis(boxsizes_vector);
            cilia.SAVAlike_CBF_measurement;
        else
            load(savename);
            cilia.load_movie;
            cilia.VariableBoxSize_Analysis(boxsizes_vector);
        end
        
        cilia.gather_results;
        save(savename, 'cilia');
        
        clearvars -global
        clear cilia savename
        
        timeelapsed(i) = toc;
        tic;
        timetogo = mean(timeelapsed(startfrom:i))*(stopat-i);
        cprintf('*[1 .3 0]',['\ntime remaning approx ',num2str(timetogo/60),' minutes\n']);
    end
catch err

	disp('Script Failed');
	disp(err.message);
	
end

disp('Script Finished')


end %function





%% ----------------------- GUI settings ---------------------- %%

function [data_folder, analysis_folder, files_string, boxsizes_vector ] =...
    GUI_settings(data_folder, analysis_folder, files_string, boxes_mode)


if isempty(data_folder)
    data_folder = pwd;
end %if

if isempty(analysis_folder)
    analysis_folder = [data_folder,'_Analysis'];
end %if

if isempty(files_string)
    files_string = '*.movie';
end %if

if isempty(boxes_mode)
    boxes_mode = 'custom';
    boxsizes_vector = [];
end

% figure
hf = figure('CloseRequestFcn',@my_close_req);
hf.Units = 'Normalized';
hf.Position = [0.25 0.1 0.5 0.8];


% panel choose data_folder
hp1 = uipanel;
hp1.Title = 'Type path to data folder or use the button to select the folder interactively';
hp1.FontSize = 10;
hp1.Units = 'normalized';
hp1.Position = [0 0.92 1 0.08];
% hp1.BorderType = 'none';

% edit choose data_folder
hedt1 = uicontrol(hp1);
hedt1.Style = 'edit';
hedt1.String = data_folder;
hedt1.FontSize = 10;
hedt1.HorizontalAlignment = 'left';
hedt1.Units = 'normalized';
hedt1.Position = [0.01 0.4 0.89 0.5];
hedt1.Callback = @type_data_folder;


% button choose data_folder
hb1 = uicontrol(hp1);
hb1.Style = 'pushbutton';
hb1.String = '...';
hb1.FontSize = 10;
hb1.HorizontalAlignment = 'center';
hb1.Units = 'normalized';
hb1.Position = [0.91 0.38 0.08 0.54];
hb1.Callback = @choose_data_folder;


% panel choose files string
hp2 = uipanel;
hp2.Title = 'Type string to identify files, press enter to update list preview';
hp2.FontSize = 10;
hp2.Units = 'normalized';
hp2.Position = [0 0.84 1 0.08];

% edit choose files string
hedt2 = uicontrol(hp2);
hedt2.Style = 'edit';
hedt2.String = files_string;
hedt2.FontSize = 10;
hedt2.HorizontalAlignment = 'left';
hedt2.Units = 'normalized';
hedt2.Position = [0.01 0.4 0.89 0.5];
hedt2.Callback = @update_list;

% panel preview files list
hp3 = uipanel;
hp3.FontSize = 10;
hp3.Units = 'normalized';
hp3.Position = [0 0.3 1 0.54];

hlist2 = uicontrol(hp3);
hlist2.Style = 'listbox';
hlist2.Units = 'normalized';
hlist2.BackgroundColor = 'w';
hlist2.HorizontalAlignment = 'left';
hlist2.String = '';
hlist2.Position = [0 0 1 1];


% panel choose analysis_folder
hp4 = uipanel;
hp4.Title = 'Type path to analysis folder or use the button to select the folder interactively';
hp4.FontSize = 10;
hp4.Units = 'normalized';
hp4.Position = [0 0.22 1 0.08];
% hp4.BorderType = 'none';

% edit choose analysis_folder
hedt4 = uicontrol(hp4);
hedt4.Style = 'edit';
hedt4.String = analysis_folder;
hedt4.FontSize = 10;
hedt4.HorizontalAlignment = 'left';
hedt4.Units = 'normalized';
hedt4.Position = [0.01 0.4 0.89 0.5];
hedt4.Callback = @update_analysis_folder;

% button choose analysis_folder
hb4 = uicontrol(hp4);
hb4.Style = 'pushbutton';
hb4.String = '...';
hb4.FontSize = 10;
hb4.HorizontalAlignment = 'center';
hb4.Units = 'normalized';
hb4.Position = [0.91 0.38 0.08 0.54];
hb4.Callback = @choose_analysis_folder;


% panel choose boxsizes
hrbg5 = uibuttongroup;
hrbg5.Title = 'Choose type of analysis';
hrbg5.FontSize = 10;
hrbg5.Units = 'normalized';
hrbg5.Position = [0 0.14 1 0.08];
hrbg5.SelectionChangedFcn = @update_boxsizes;

% hp5.BorderType = 'none';

% radio button group
hrb5(1) = uicontrol(hrbg5);
hrb5(1).Style = 'radiobutton';
hrb5(1).String = 'quick';
hrb5(1).FontSize = 10;
hrb5(1).Units = 'Normalized';
hrb5(1).Position = [0.01 0.35 0.1 0.5];

hrb5(2) = uicontrol(hrbg5);
hrb5(2).Style = 'radiobutton';
hrb5(2).String = 'sigmoids';
hrb5(2).FontSize = 10;
hrb5(2).Units = 'Normalized';
hrb5(2).Position = [0.11 0.35 0.1 0.5];

hrb5(3) = uicontrol(hrbg5);
hrb5(3).Style = 'radiobutton';
hrb5(3).String = 'custom';
hrb5(3).FontSize = 10;
hrb5(3).Units = 'Normalized';
hrb5(3).Position = [0.24 0.35 0.1 0.5];

hedt5 = uicontrol(hrbg5);
hedt5.Style = 'edit';
hedt5.FontSize = 10;
hedt5.HorizontalAlignment = 'left';
hedt5.Units = 'normalized';
hedt5.Position = [0.35 0.35 0.64 0.5];
hedt5.Callback = @type_custom_boxsizes;

% start analysis
hpb6 = uicontrol(hf);
hpb6.Style = 'pushbutton';
hpb6.String = 'Start the analysis';
hpb6.FontSize = 14;
hpb6.FontWeight = 'b';
hpb6.BackgroundColor = 'g';
hpb6.Units = 'normalized';
hpb6.Position = [0.4 0.02 0.2 0.1];
hpb6.Callback = @start_analysis;

% initialise
hrbg5.SelectedObject = ...
    hrb5( arrayfun(@(i) strcmp(boxes_mode, hrb5(i).String), 1:numel(hrb5) ) );
type_data_folder;
update_boxsizes;
abort = false;

% output
waitfor(hf)
if abort
    error('User close the settings window - execution stopped.');
end


    function [] = type_data_folder(~,~)
        % check that the folder exists, updates the automatic analysis
        % folder name and preview list
        
        % check the folder exists
        if ~isdir(hedt1.String)
            hedt1.BackgroundColor = 'r';
            return
        else
            hedt1.BackgroundColor = 'w';
        end %if
        
        data_folder = hedt1.String;
        
        if strcmp(data_folder(end),'/') || strcmp(data_folder(end),'\')
            data_folder = data_folder(1:end-1);
            hedt1.String = data_folder;
        end
        
        % update analysis folder name
        hedt4.String = [data_folder,'_Analysis'];
        
        update_list;
        update_analysis_folder;
        
    end %function


    function [] = choose_data_folder(~,~)
        % choose with interface where to take daata from
        
        data_folder = uigetdir(data_folder);
        
        % update editable field back
        hedt1.String = data_folder;
        
        type_data_folder;
        
    end %function


    function [] = update_list(~,~)
        
        % retrieve string
        tempstr = hedt2.String;
        data_folder = hedt1.String;
        
        % create list
        fl = dir(fullfile(data_folder,tempstr));
        arraylist = arrayfun(@(i)fl(i).name,(1:numel(fl))','UniformOutput',0);
        
        hlist2.String = arraylist;
        files_string = hedt2.String;
    end %function


    function [] = update_analysis_folder(~,~)
        
        % two cases, either it exists or need to be created
        
        switch isdir(hedt4.String)
            
            case false
                status = mkdir(hedt4.String);
                rmdir(hedt4.String);
                
                if ~status
                    hedt4.BackgroundColor = 'r';
                else
                    hedt4.BackgroundColor = 'g';
                end %if
                
                
            case true
                
                alphabet = char(uint8([48:57,97:122]));
                random_filename = fullfile(hedt4.String,[alphabet(randi(numel(alphabet),1,20)),'.txt']);
                
                [fid,errmsg] = fopen(random_filename, 'w');
                if ~isempty(errmsg)
                    hedt4.BackgroundColor = 'r';
                else
                    hedt4.BackgroundColor = 'g';
                    fclose(fid);
                    delete(random_filename);
                end %if
                
        end %switch
        
        analysis_folder = hedt4.String;
        
        
    end %function

    function [] = choose_analysis_folder(~,~)
        % choose with interface where to take daata from
        
        analysis_folder = uigetdir(analysis_folder);
        
        % update editable field back
        hedt4.String = analysis_folder;
        
        update_analysis_folder;
        
    end %function


    function [] = update_boxsizes(~,~)
        
        
        boxes_mode = hrbg5.SelectedObject.String;
        
        switch hrbg5.SelectedObject.String
            
            case 'quick'
                boxsizes_vector = get_boxsizes_vector(boxes_mode);
                hedt5.Enable = 'off';
            case 'sigmoids'
                boxsizes_vector = get_boxsizes_vector(boxes_mode);
                hedt5.Enable = 'off';
            case 'custom'
                hedt5.Enable = 'on';
                if isempty(hedt5.String)
                    boxsizes_vector = get_boxsizes_vector(boxes_mode);
                end %if
        end %switch
        
        hedt5.String = num2str(boxsizes_vector,'%d  ');
        
    end %function


    function [] = type_custom_boxsizes(~,~)
        % this can only be done if radiobutton is on custom mode
        
        boxsizes_vector = str2num(hedt5.String); %#ok<ST2NM>
        
    end %function

    function [] = start_analysis(~,~)
        
        if ~isdir(analysis_folder)
            mkdir(analysis_folder)
        end %if
        
        delete(hf)
        
    end %function

    function [] = my_close_req(~,~)
        
        abort = true;
        delete(hf)
        
    end %function

end %function

%% ----------------------- Additional functions ---------------------- %%


function [boxsizes_vector] = get_boxsizes_vector(boxes_mode)

boxsizes_vector_quick = [32 64 128 256 512 1024];
boxsizes_vector_sigmoids = [16 32 48 64 96 128 160 192 224 256 340 512 1024];

switch boxes_mode
    case 'quick'
        boxsizes_vector = boxsizes_vector_quick;
    case 'sigmoids'
        boxsizes_vector = boxsizes_vector_sigmoids;
    case 'custom'
        boxsizes_vector = boxsizes_vector_quick;
    otherwise
        error('Unknown analysis mode.');
end

end %function




function [] = splash()
%%
home


% choose characters

ULc = ' ';   % upper left corner
URc = ' ';   % upper right corner
LLc = '|';   % lower left corner
LRc = '|';   % lower right corner
HL  = '_';   % horizontal line
VL  = '|';   % vertical line
 
% print

fprintf('\n\n');
fprintf('%s', ULc);
fprintf('%s', repmat(HL, 1, 69));
fprintf('%s', URc);
fprintf('\n');

fprintf('%s', VL);
fprintf('%s', repmat(' ',1,69));
fprintf('%s', VL);
fprintf('\n');

fprintf('%s', VL);
fprintf('%s', repmat(' ',1,19));
fprintf('%s', 'Semi-Automatic ALI DDM Analysis');
fprintf('%s', repmat(' ',1,19));
fprintf('%s', VL);
fprintf('\n');

fprintf('%s', VL);
fprintf('%s', repmat(' ',1,24));
fprintf('%s', 'Based on DDM_Analysis');
fprintf('%s', repmat(' ',1,24));
fprintf('%s', VL);
fprintf('\n');

fprintf('%s', VL);
fprintf('%s', repmat(' ',1,69));
fprintf('%s', VL);
fprintf('\n');

fprintf('%s', VL);
fprintf('%s', repmat(' ',1,5));
fprintf('%s', [char(169),' Luigi Feriani et al., 2017']);
fprintf('%s', repmat(' ',1,36));
fprintf('%s', VL);
fprintf('\n');

fprintf('%s', VL);
fprintf('%s', repmat(' ',1,69));
fprintf('%s', VL);
fprintf('\n');

fprintf('%s', VL);
fprintf('%s', repmat(' ',1,6));
fprintf('%s', 'This work is released under license CC BY-NC-ND 4.0');
fprintf('%s', repmat(' ',1,12));
fprintf('%s', VL);
fprintf('\n');

fprintf('%s', VL);
fprintf('%s', repmat(' ',1,69));
fprintf('%s', VL);
fprintf('\n');

fprintf('%s', LLc);
fprintf('%s', repmat(HL, 1,69));
fprintf('%s', LRc);
fprintf('\n\n');

end %function










