addpath('./Tworit_functions_toolbox/');
%% Deliverable#4 Tworit v1.0 S Plotting Geometry from saved Geometry in .mat files

% Inputs appear as the code is run:

%% ==================================================== Plot Section =======================================================

Plot_enable = input('Do you want to plot the Geometry from a mat file? : 1 for Yes, 0 for No: [0]: ');

if isempty(Plot_enable) || (Plot_enable ~= 1)
    Plot_enable = 0;
end

if Plot_enable

disp('===========================================================');
disp('Plot section: ');

disp('========================File Selection===================================');


file_selection = 1;

while file_selection

[file_name,selected_dir] = uigetfile('*.mat');



if isequal(file_name,0)
   disp('User selected No file: ');
   disp('Select file again: ')
else
   disp(['User selected ', fullfile(selected_dir,file_name)]);
   file_selection = 0;
end

end

plot_info.file = [selected_dir, '/', file_name];

GEO = load(plot_info.file);


Plot_GEO(GEO.geo);

end
