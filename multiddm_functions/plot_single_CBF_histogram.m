function [ hf ] = plot_single_CBF_histogram( MergedData, BoxSizeIndex, binedges, flag_numbers )
%plot_single_CBF_histogram Takes the structure MergedData
%and plots the CBF in a histogram.
%   MergedData should be a scalar structure from merge_SampleType_data
%   or it can be an array of such structure. I'm assuming that the
%   sampletype does not change, and that the only thigs that changes
%   between entries of the structure array is the timepoints
%   
%   BoxSizeIndex should be the index corresponding to the tile size of
%   which you want to plot the histogram. Defaults to finding the tile with
%   size closest to 10um x 10um
%   
%   binedges is an array with the edges of the bins for the histogram
%   (similar to what MATLAB's histogram function accept as input). Defaults
%   to plotting a histogram with bins 1Hz wide, between 0 and 30
%
%   flag_numbers: if True, wites mean +- std in the legend
% 


%{
% Version 1.0
% © Luigi Feriani 2019 (luigi.feriani@gmail.com) 
% 
% plot_single_CBF_histogram.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% 
% Original work:
% 
% Chioccioli, M.*, Feriani, L.*, Kotar, J., Bratcher, P. E.**, Cicuta, P.**, Nature Communications 2019
% "Phenotyping ciliary dynamics and coordination in response to CFTR-modulators 
% in Cystic Fibrosis respiratory epithelial cells"
%}

%% input check

% first check that the conditions for MergedData have been met

% check that MergedData is a scalar
if ~isscalar(MergedData)
    error('MergedData can only be a struct, not an array of structs.')
end

% now the rest of the input check

% flag numbers
if nargin < 4 || isempty(flag_numbers)
    flag_numbers = false;
end


% binedges
if nargin < 3 || isempty(binedges) || ~isvector(binedges)
    binedges = 0:1:30;
end %if


% boxsize index
if nargin < 2 || isempty(BoxSizeIndex) || BoxSizeIndex > numel(MergedData(1).window_area_um2)
    [~,bsi] = min(abs(100 - MergedData.window_area_um2)); % find the closest to 100um2
else
    bsi = BoxSizeIndex;
end %if



%% automatic process fix parameters

binedges = sort(binedges);
binedges(binedges < 0 | binedges > 60) = [];

cmap = lines(1);

%% plot

% prepare figure
hf = figure;

% prepare axes
ha = axes;
% ha.Position = [0.15 0.15 0.8 0.75];
ha.XLim = minmax(binedges);
ha.Box = 'on';
ha.NextPlot = 'add';

% labels
ha.XLabel.String = 'CBF, [Hz]';
ha.XLabel.FontSize = 16;
ha.YLabel.String = 'Counts, normalized';
ha.YLabel.FontSize = 16;
ha.Title.String = [sprintf('%s    ',MergedData.sampletype_str),...
    sprintf('%s    ',MergedData.timepoint_str),...
    sprintf('%s    ',MergedData.donors_str{:}),...
    sprintf('%s    ',MergedData.inserts_str{:}),...
    sprintf('%s    ',MergedData.positions_str{:})];ha.Title.Interpreter = 'none';
ha.Title.FontSize = 16;

% actual plot

% histogram with edge
hhs = histogram(ha, MergedData.Frequency_Hz{bsi}, binedges);
hhs.Normalization = 'Probability';
hhs.DisplayStyle = 'Stairs';
hhs.LineWidth = 1.5;
hhs.EdgeColor = cmap;

% edgeless histogram
hhb = histogram(ha, MergedData.Frequency_Hz{bsi}, binedges);
hhb.Normalization = 'Probability';
hhb.EdgeColor = 'none';
hhb.FaceColor = cmap;
hhb.FaceAlpha = 0.2;

% legend
if flag_numbers
    leg = ['CBF = (',...
        num2str(nanmean(MergedData.Frequency_Hz{bsi}),3),...
        ' ± ',num2str(nanstd(MergedData.Frequency_Hz{bsi}),3), ') Hz'];
end %if



% legend
hleg = legend(ha, hhb, leg);
hleg.Box = 'off';
hleg.FontSize = 10;
hleg.Interpreter = 'none';



end

