function [] = plot_single_sigmoid_errorbar(MergedData, flag_confint, whichshoulder)
%plot_single_sigmoid_errorbar takes a MergedData structure (as output from
%merge_SampleType_data) and plots a single sigmoid with the fitting line
%and errorbar

%{
% Version 1.0
% © Luigi Feriani 2019 (lf352@cam.ac.uk, luigi.feriani@gmail.com) 
% 
% plot_single_sigmoid_errorbar.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.
% 
% Original work:
% 
% Chioccioli, M.*, Feriani, L.*, Kotar, J., Bratcher, P. E.**, Cicuta, P.**, Nature Communications 2019
% "Phenotyping ciliary dynamics and coordination in response to CFTR-modulators 
% in Cystic Fibrosis respiratory epithelial cells"
%}

if nargin < 3 || isempty(whichshoulder) || ~any(strcmp(whichshoulder,{'left','right','both'}))
    whichshoulder = 'right';
end 
    

if nargin < 2 || isempty(flag_confint)
    flag_confint = false;
end

% prepare figure
hf = figure;

% prepare axes
ha = axes;
box on;
hold on;

% axes limits
setsemilogx
ha.YLim = [0 max(vertcat(MergedData.Damping_Hz{:}))];
ha.XLim = [1e0 1e5];
ha.XTick = logspace(0,5,6);

% errorbar
he = errorbar(MergedData.window_area_um2, MergedData.med_Damping_Hz,...
    MergedData.ler_Damping_Hz, MergedData.uer_Damping_Hz);
he.LineWidth = 1.2;
he.LineStyle = 'none';
he.Marker = 'o';
he.MarkerFaceColor = 'w';

% fit
xx = logspace(0,5,1e3);
hpf = plot(xx, MergedData.Damping_Hz_fit_out2(log10(xx)));
hpf.LineWidth = 1.2;
hpf.Color = ha.ColorOrder(1,:);
hpf.Tag = 'fitline';

% vertical line - is now the shoulder
switch whichshoulder
    case 'right'
        hpv = plot(10.^(MergedData.Damping_Hz_fit_out2.mu)*[1 1].*exp(2), [0 20*ha.YLim(2)]);
        hpv.Tag = 'rightline';
    case 'left'
        hpv = plot(10^(MergedData.Damping_Hz_fit_out2.mu)*[1 1]./exp(2), [0 20*ha.YLim(2)]);
        hpv.Tag = 'leftline';
    case 'both'
        hpv(1) = plot(10.^(MergedData.Damping_Hz_fit_out2.mu)*[1 1].*exp(2), [0 20*ha.YLim(2)]);
        hpv(1).Tag = 'rightline';        
        hpv(2) = plot(10^(MergedData.Damping_Hz_fit_out2.mu)*[1 1]./exp(2), [0 20*ha.YLim(2)]);
        hpv(2).Tag = 'leftline';
    otherwise
end
set(hpv,'Color', hpf.Color);
set(hpv,'LineStyle', hpf.LineStyle);
set(hpv,'LineWidth', hpf.LineWidth);


if flag_confint
    dummy = 10.^(par_confint(MergedData.Damping_Hz_fit_out2,'mu',0.68));
    
    switch whichshoulder
        case 'right'
            hpp = patch('XData',dummy([1 2 2 1]).*exp(2),'YData', ha.YLim([1 1 2 2]));
            hpp.Tag = 'rightpatch';
        case 'left'
            hpp = patch('XData',dummy([1 2 2 1])./exp(2),'YData', ha.YLim([1 1 2 2]));
            hpp.Tag = 'leftpatch';
        case 'both'
            hpp(1) = patch('XData',dummy([1 2 2 1]).*exp(2),'YData', ha.YLim([1 1 2 2]));
            hpp(1).Tag = 'rightpatch';
            hpp(2) = patch('XData',dummy([1 2 2 1])./exp(2),'YData', ha.YLim([1 1 2 2]));
            hpp(2).Tag = 'leftpatch';
        otherwise
    end
    set(hpp,'FaceColor', hpf.Color);
    set(hpp,'FaceAlpha', 0.2);
    set(hpp,'EdgeColor', 'none');
end

 
% labels and titles
ha.XLabel.String = 'DDM Window Area, [\mum^2]';
ha.YLabel.String = '\tau_c^{-1}, [s^{-1}]';
ha.XLabel.FontSize = 16;
ha.YLabel.FontSize = 16;
ha.Title.String = [sprintf('%s    ',MergedData.sampletype_str),...
    sprintf('%s    ',MergedData.timepoint_str),...
    sprintf('%s    ',MergedData.donors_str{:}),...
    sprintf('%s    ',MergedData.inserts_str{:}),...
    sprintf('%s    ',MergedData.positions_str{:})];
ha.Title.Interpreter = 'none';
ha.Title.FontSize = 16;

end %function

%%%%%%%%%%%%%%

        
        