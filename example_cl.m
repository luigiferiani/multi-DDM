cd %insert_data_folder_path


clear all
close all
home

filelist = dir('40X_BF*.movie'); %change this as needed

savepath = %insert_analysis_folder;
mkdir(savepath)

tic;
timeelapsed = zeros(numel(filelist),1);
startfrom = 1;
stopat = numel(filelist);

try
    for i = startfrom:stopat 
        
        cprintf('*[0 0 .4]',['\n',filelist(i).name,'\n'])
		
		
        savename = filelist(i).name;
        savename = savename(1:end-6);  % modify this to suit other files extensions
        savename = fullfile(savepath,[savename,'.mat']);
        
        % load(savename);
        
        cilia = DDM_Analysis(filelist(i).name);
        cilia.N_couple_frames_to_average = 200;
%         cilia.VariableBoxSize_Analysis([32 64 128 256 512 1024]); 
        cilia.VariableBoxSize_Analysis([16 32 48 64 96 128 160 192 224 256 340 512 1024]); 
        cilia.SAVAlike_CBF_measurement;
        
        %         cilia.load_movie;
        %         cilia.VariableBoxSize_Analysis([32 64]);
        
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

	disp('Script failed')
	err.message
end



