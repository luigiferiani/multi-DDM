% Example script to plot all data aggregated by populate_DDM_Sigmoids_struct

% clear workspace
cch

% load aggregated data
aggregated_data = 'Exp3_data.mat';
load(aggregated_data,...
    'SampleType',...
    'figures_folder',...
    'sampletypes',...
    'timepoints',...
    'donors',...
    'inserts',...
    'positions',...
    'boxsizes_vector',...
    'q_limits_1oum');

%% Simple use of the MergedData structure:

% plot a sigmoidal curve and a CBF distribution for data from each
% sampletype, timepoint, insert, position, donors combination

% in this case, we have different donors, and inserts, that we want to keep
% separate, but we do not have positions that we want to keep separate.
        
for stc = 1:numel(sampletypes) % sampletype counter
    for tpc = 1:numel(timepoints) % timepoint counter
        for dc = 1:numel(donors) % donor counter
            for ic = 1:numel(inserts) % insert counter
                
                % output stuff
                disp_C = {sampletypes{stc},...
                    timepoints{tpc},...
                    donors{dc},...
                    inserts{ic}};
                disp_str = strjoin(disp_C,', ');
                disp('Plotting data from: ')
                disp(disp_str)
                
                % invoke function that extracts data from the SampleType
                % structure for easier plotting
                MergedData = merge_SampleType_data(SampleType,...
                    stc, tpc, donors{dc}, inserts{ic}, '');
                
                % plot sigmoidal curve with left shoulder and shading
                % showing the confidence interval
                plot_single_sigmoid_errorbar(MergedData, true, 'left')
                
                % plot CBF distribution measured at the tile size closest
                % to 10um x 10um (default), using 0.5Hz wide bins between 0
                % and 30Hz, writing mean and std in the legend
                plot_single_CBF_histogram(MergedData, [], 0:0.5:30,1)

%                 structstruct(MergedData)

            end
        end % for dc
    end % for tpc
end %for stc


% It is also possible to use merge_SampleType_data to pool together data
% taken on different inserts, or at different positions, or even from
% different donors. Just pass an array as the relevant input:
% if one wanted to pool together all donors and inserts:
% MergedData =  merge_SampleType_data(SampleType, stc, tpc, donors{dc}, inserts, '');
% if one wanted to pool together all donors and inserts:
% MergedData =  merge_SampleType_data(SampleType, stc, tpc, donors, inserts, '');
% (of course, this requires tweaking the for loops above)
        
% It is then possible to use an array of MergedData structures to create
% figures where multiple sigmoids or multiple CBF distributions are
% overlaid, but that means writing code that is quite specialised to the
% task/experiment, and is beyond the scope of this commit
