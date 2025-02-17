addpath('./Tworit_functions_toolbox/');
%% Deliverable#2 Tworit v1.0 S Plotting S Parameters from saved S parameters in .mat files

% Inputs appear as the code is run:

%% ==================================================== Plot Section =======================================================

Plot_enable = input('Do you want to plot the S parameter you like now? : 1 for Yes, 0 for No: [0]: ');

if isempty(Plot_enable) || (Plot_enable ~= 1)
    Plot_enable = 0;
end

plot_info.same_diff = 0;

while Plot_enable

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

MM = load(plot_info.file);


Nbn = size(MM.b.SRR, 2);
Ntn = size(MM.b.STT, 2);



disp("T is the top waveguide: R is the base Waveguide: ")

disp("Select the S parameter you want: ")

S_valid = 1;

while S_valid

  Config_S = input('1 for SRR, 2 for STT, 3 for SRT, 4 for STR: ');
  
  if (Config_S == 1) || (Config_S == 2) || (Config_S == 3) || (Config_S == 4) 
      S_valid  = 0;
  end
  
  if S_valid == 1
      disp('Given configuration is wrong: Please input one of the following values only: ');
  end
  
end

mode_conf = 1;

while mode_conf

 
if Config_S == 1
  disp('Check the serial number of the modes in MODES.txt before entering mode serial numbers below');
      
  Config_m = input('Serial number of first mode of R: ');
  Config_n = input('Serial number of second mode of R: ');
  
  if (Config_m < 1) || (Config_n < 1) || (Config_m > Nbn) || (Config_n > Nbn)
     mode_conf = 1;  
  else
     mode_conf = 0;
  end
     
elseif Config_S == 2
    
  disp('Check the serial number of the modes in MODES.txt before entering mode serial numbers below');
      
  Config_m = input('Serial number of first mode of T: ');
  Config_n = input('Serial number of second mode of T: ');
  
  if (Config_m < 1) || (Config_n < 1) || (Config_m > Ntn) || (Config_n > Ntn)
     mode_conf = 1;  
  else
     mode_conf = 0;
  end
      
elseif Config_S == 3
   disp('Check the serial number of the modes in MODES.txt before entering mode serial numbers below');
      
  Config_m = input('Serial number of first mode of R: ');
  Config_n = input('Serial number of second mode of T: ');
  
  if (Config_m < 1) || (Config_n < 1) || (Config_m > Nbn) || (Config_n > Ntn)
     mode_conf = 1;  
  else
     mode_conf = 0;
  end
else
   disp('Check the serial number of the modes in MODES.csv before entering mode serial numbers below');
      
  Config_m = input('Serial number of first mode of T: ');
  Config_n = input('Serial number of second mode of R: ');
  
  if (Config_m < 1) || (Config_n < 1) || (Config_m > Ntn) || (Config_n > Nbn)
     mode_conf = 1;  
  else
     mode_conf = 0;
  end
end
  
if mode_conf == 1
    disp('Configuration of modes is not valid, try again');
    disp('Possible things which are wrong can be: ');
    disp('Either the mode serial number(s) are less than 1: OR ');
    disp('The serial number(s) are out of bounds: ');
    disp(['For base waveguide (R): mode number has to be between 1 and ', num2str(Nbn), ' :']);
    disp(['For top waveguide (T): mode number has to be between 1 and ', num2str(Ntn), ' :']);
    disp('Try to enter the serial number again: ');
end

end

plot_info.Config_S = Config_S;
plot_info.Config_m = Config_m;
plot_info.Config_n = Config_n;
% plot_info.f_axis = input_info.MODES.f_axis;
% plot_info.f_axis = MM.b.f_axis;

Plot_S(plot_info);

Plot_enable = input('Want to plot a new S parameter?: 1 for Yes and 0 for No: [0]: ');
if isempty(Plot_enable) || (Plot_enable ~= 1)
 Plot_enable = 0;
end

if Plot_enable == 1
    same_diff = input('want on the same plot or a different plot? 1 for same, 0 for different: [0]: ');
    if isempty(same_diff) || (same_diff ~= 1)
        same_diff = 0;
    end
    
    plot_info.same_diff = same_diff;

end

end