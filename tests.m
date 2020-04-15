clear all
close all
boxsizes_vector = [64 512 1024];
boxsizes_encore = [128];

%% test OME
filename = 'Health Cells.czi';
filepath = 'F:\Data\Cilia_NSW\example_czi';

cilia = DDM_Analysis(filename, filepath);
cilia.set_temperature(1);
cilia.N_couple_frames_to_average = 200;
cilia.VariableBoxSize_Analysis(boxsizes_vector);
cilia.SAVAlike_CBF_measurement;
cilia.gather_results;

clearvars -global fs
cilia.VariableBoxSize_Analysis(boxsizes_encore);
cilia.gather_results;


%% test movie
clear global
clearvars cilia

filename = '40X_BF_Nm_d3.12Apr2017_15.52.18.movie';
filepath = 'D:\Data\Cilia\Data\2017_04_12';

cilia = DDM_Analysis(filename, filepath);
cilia.set_temperature(1);
cilia.N_couple_frames_to_average = 100;
cilia.VariableBoxSize_Analysis(boxsizes_vector);
cilia.SAVAlike_CBF_measurement;
cilia.gather_results;

clearvars -global fs
cilia.VariableBoxSize_Analysis(boxsizes_encore);
cilia.gather_results;

%% test avi
clear global
clearvars cilia

filename = '18.avi';
filepath = 'Z:\public_html\out\Cilia_example_videos';

cilia = DDM_Analysis(filename, filepath);
cilia.set_temperature(1);
cilia.N_couple_frames_to_average = 100;
cilia.VariableBoxSize_Analysis(boxsizes_vector);
cilia.SAVAlike_CBF_measurement;
cilia.gather_results;

clearvars -global fs
cilia.VariableBoxSize_Analysis(boxsizes_encore);
cilia.gather_results;
