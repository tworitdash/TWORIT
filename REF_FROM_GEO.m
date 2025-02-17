%% Deliverable#8 Tworit v1.0 S Parameters from a given Geometry mat file
%
% THIS SCRIPT IS FOR REFLECTION COEFFICIENTS OF CONICAL HORN ANTENNA
%
% DESCRIPTION: 
% To find out reflection coefficients of a cascaded
% cylindrical waveguide structure/ conical horn structure


% SETTINGS:
% Setting ¨a¨ is all about the geometry of the waveguide/ horn and setting
% ¨b¨ is all about the waveguide modes
% ¨c¨ is return loss based on a given mode of excitation at the base
% element of the structure 

%        a) GEOMETRY: Geometry is given from a mat file: A geometry can be
%        created by GEO_Con.m file
%           

% OUTPUT: 
%        a) Reflection/Transmission Coefficients of mode pairs in 3D
%        format:It can be saved with an option

%           DIM1: Frequency axis
%           DIM2: Modes at the base waveguide
%           DIM3: Modes at the top waveguide
%% 
clear;
close all;

addpath('./Tworit_functions_toolbox/');
change_options = 1;
while change_options
  %% ==================================================== Geometry of the Antenna ============================================================
    disp('===========================================================');
    disp('Geometry section: ');

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

    geo_info.file = [selected_dir, '/', file_name];

    GEO = load(geo_info.file);
    input_info.GEO = GEO.geo;
    
    input_info.Plot_Geometry = input('Do you want to plot the Geometry?: Select 1 for YES, Select 2 for NO, Default is NO: ');
    
    if isempty(input_info.Plot_Geometry)
        input_info.Plot_Geometry = 0;
    end
    
    if input_info.Plot_Geometry == 1
        Plot_GEO(input_info.GEO);
    end
    %% ========================================== Frequency and modes ==============================================
   f_invalid = 1;
    
    while f_invalid
        input_info.MODES.f_begin = input('Input the start frequency of the operation in [GHz]: ');
        input_info.MODES.f_end = input('Input the end frequency of the operation in [GHz]: ');
        
        if input_info.MODES.f_begin < input_info.MODES.f_end
            f_invalid = 0;
        end
        
        if f_invalid == 1
            disp('Invalid choice of frequencies, make sure the end frequency is greater than the start frequency: ');
        end
       
       
    end
    
    input_info.MODES.Nf = input('Input the number of frequency points you need: ');
    
     
    f_cutoff = fc(input_info.GEO.E(1).r, 1, 1);
    f_cutofft = fc(input_info.GEO.E(end).r, 1, 1);
        
    if f_cutoff(1) > (input_info.MODES.f_begin * 1e9)
        disp('The tool will discard the frequencies from : '); disp(input_info.MODES.f_begin); disp(' GHz till '); disp(f_cutoff(1) * 1e-9); disp('GHz'); 
        disp(' Because they are evanascent for the first waveguide in the Geometry: ');
        input_info.MODES.f_axis = linspace(f_cutoff(1), input_info.MODES.f_end * 1e9, input_info.MODES.Nf);
    else
        input_info.MODES.f_axis = linspace(input_info.MODES.f_begin * 1e9, input_info.MODES.f_end * 1e9, input_info.MODES.Nf);
    end
    
    disp(' Smart mode selection strategy starts: ');
    
   
        modes = load('Xmn.mat');
        modes = modes.Xmn;
       
        [~, Nb] = find(f_cutoff(input_info.MODES.f_axis(end) > f_cutoff));
        [~, Nt] = find(f_cutoff(input_info.MODES.f_axis(end) > f_cutofft));
        
        Nbn = length(Nb);
        Ntn = length(Nt);
        input_info.MODES.Nbn = Nbn;
        input_info.MODES.Ntn = Ntn;
        
        if round(Ntn/Nbn) < 25
            Nbn = 20;
            Ntn = round(Ntn/Nbn) .* Nbn;
            
        end
         disp(' Smart mode selection strategy ends: ');
        
         print_modes = input('Want to print the modes? 1 for YES, 0 for NO [0]: ');
        
        if isempty(print_modes)
            print_modes = 0;
        end
        
        if print_modes == 1
            disp('====================================================');
            disp('PORT 1: TOP WAVEGUIDE');
            for m = 1:Nbn
               disp(modes(m).mode); disp(modes(m).m); disp(modes(m).n); disp('Orientation: '); disp(modes(m).pol * 180/pi); disp(' [deg]'); 
            end
            



            disp('====================================================');
            disp('PORT 2: BASE WAVEGUIDE');
            
            for m = 1:Ntn
                disp(modes(m).mode); disp(modes(m).m); disp(modes(m).n); disp('Orientation: '); disp(modes(m).pol * 180/pi); disp(' [deg]'); 
            end
            
        end 

        
     
    
    
    
    change_options = input('Give 0 is you agree with this options, 1 if you want to reintroduce them [0]:');
    if isempty(change_options)
        change_options = 0;
    end
end

%% =========================== Calculation of S Parameters ===========================================================


disp('Calling the S paremeter function: ');




[RLRR_TE11, Srr, Stt, Str, Srt] = GSM_v1(input_info);




disp('S Parameters Computed: ');


% figure; plot(input_info.MODES.f_axis .* 1e-9, db(STT(:, 1, 1)));

b.STT = Stt;
b.SRR = Srr;
b.STR = Str;
b.SRT = Srt;
b.RLRR_TE11 = RLRR_TE11;
b.f_axis = input_info.MODES.f_axis;

%% ============================================== Save Data ==============================================================

%% S Parameters

disp("Save Data Section: ")

save_ = input('Do you want to save the results in a file? 1 for Yes, 0 for No: [1]: ');

if isempty(save_)
    save_ = 0;
end



if (save_ == 1)
    disp('Selct the folder in the UI: ')
    selected_dir = uigetdir();
    file_name = input("Input the file name you want: Use single quotes when writing file name such as: 'Test': ");
   
    save([selected_dir, '/', file_name], 'b');
end


disp(['Data saved in', fullfile(selected_dir,file_name), '.mat']);


%% end ===================================================================================================================
