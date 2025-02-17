%% Deliverable#6 Tworit v1.0 Far Fields From saved Geometry and General Scattering Matrix [S parameters]
%
% THIS SCRIPT IS FOR FAR FIELDS OF CONICAL HORN ANTENNA WHEN THE GEOMETRY
% AND THE S PARAMETER MATRIX OR THE GENERAL SCATTERING MATRIX ARE TAKEN
% FROM SAVED FILES
%
% DESCRIPTION: 
% To find out far fields of a cascaded
% cylindrical waveguide structure/ conical horn structure


% SETTINGS:
% Setting ¨a¨ is all about the geometry of the waveguide/ horn and setting
% ¨b¨ is all about the waveguide modes
% ¨c¨ is all about far fields inputs

%        a) GEOMETRY: Geometry is a input here from saved mat files. 

%        b) S PARAMETERS/ GENERAL SCATTERING MATRIX: GENERAL SCATTERING
%        MATRIX IS CHOSEN FROM A SAVED MAT FILE.

% OUTPUT: 
%       
%        b) Far Fields
%           DIM1: Frequency axis
%           DIM2: Azimuth
%           DIM3: Elevation
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
    
    disp('===========================================================');
    disp('GSM section: ');

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

    gsm_info.file = [selected_dir, '/', file_name];

    MM = load(gsm_info.file);
    
    input_info.FF.GSM = MM.b;
    
    Nbn = size(MM.b.SRR, 1);
    input_info.MODES.Nbn = Nbn;
    
    Ntn = size(MM.b.STT, 1);
    input_info.MODES.Ntn = Ntn;
    
    modes = load('Xmn.mat');
    modes = modes.Xmn;
    f_cutoff = fc(input_info.GEO.E(1).r, 1, 1);
    
    mff_in = 1;
       while (mff_in)
            disp('Far field Options:')
            
            mff = input(['Which modes excitation you want at the base waveguide: Input something from 1 to ', num2str(Nbn), ', Ex: [1 2 3] , Default is [1]: Fundamental mode: ']);
            
            if isempty(mff)
               mff = 1;
            end
            
            input_info.FF.ar = input(['Relative excitation amplitudes, ex: [1 0 0], has to be the same dimension as the excitation modes input given before: [1] ']);
            
            if isempty(input_info.FF.ar)
                input_info.FF.ar = ones(length(mff));
            end
            
            if (min(mff) >= 1) && (max(mff) <= Nbn) && (length(mff) == length(input_info.FF.ar))
                mff_in = 0;
                disp(['The excitation chosen is/are the mode(s) ']);
                for mffi = 1:length(mff)
                    disp(modes(mffi).mode); disp(modes(mffi).m); disp(modes(mffi).n); disp('Orientation: '); disp(modes(mffi).pol * 180/pi); disp(' [deg]'); 
                end
                input_info.FF.ex_mode = modes(mff);
                input_info.FF.mff = mff;
                
                % Frequency Selection =====================================
                disp(['Frequencies should be in between ', num2str(f_cutoff(min(mff))*1e-9), 'GHz and ', num2str(input_info.FF.GSM.f_axis(end)*1e-9), 'GHz']);
                
                fopt = 1;
                while fopt 
                
                    input_info.FF.fopt = input('Do you want to consider all the frequency points in this range? 1 YES, 0 No [1]:');
                    if isempty(input_info.FF.fopt)
                        input_info.FF.fopt = 1;
                    end
                    
                    if (input_info.FF.fopt == 0) || (input_info.FF.fopt == 1)
                        fopt = 0;
                        [~, fidx] = min(abs(input_info.FF.GSM.f_axis - f_cutoff(min(mff))));
                        if input_info.FF.GSM.f_axis(fidx) < f_cutoff(min(mff))
                            fidx = fidx + 1;
                        end
                        
                        f_axis_ff = input_info.FF.GSM.f_axis(fidx:end);
                        if (input_info.FF.fopt == 1)
                            input_info.FF.Nfi = 1:length(f_axis_ff);
                            input_info.FF.f_axis = f_axis_ff;
                        else
                            disp(['Manually select the frequency indices: ']);
                            disp(['Frequencies selected for S parameters are in (Hz): ']);
                
                            for fi = 1:length(f_axis_ff)
                                disp('Serial Number (index): '); disp(fi); disp('Frequency in Hz: '); disp(f_axis_ff(fi)); 
                            end
                
                            Nfi = 1;
                
                            while Nfi
                                input_info.FF.Nfi = input('input the indices you want: ex: [1 2 5 8]: default [1]');
                                if isempty(input_info.FF.Nfi)
                                    input_info.FF.Nfi = 1;
                                end
                                if min(input_info.FF.Nfi >= 1) && min(input_info.FF.Nfi <= length(f_axis_ff))
                                    Nfi = 0;
                                    input_info.FF.f_axis = f_axis_ff(input_info.FF.Nfi);
                                else
                                    disp('Indices are wrong, please type in again: ');
                                    Nfi = 1;
                                end
                            end
                        end
                        
                    else
                        disp(['This option is wrong, try inputing again, should be 0 or 1: ']);
                        fopt = 1;
                    end
                end
                
               
                
                
                % Normalization selection =================================
                
                norm = 1;
                while norm
                    
                    input_info.FF.norm = input('Normalization: 1 for Power normalization, 2 for FEKO style mormalization: [2] ');
                        if isempty(input_info.FF.norm)
                            input_info.FF.norm = 2;
                        end
                    
                    if (input_info.FF.norm == 1) || (input_info.FF.norm == 2)
                        norm = 0;
                    else
                        disp(['Invalid normalization input, try to input again: ']);
                        norm = 1;
                    end
                end
               
                % Aperture reflection selection ===========================
                
                ap = 1;
                while ap
                    disp(['This option is for open aperture, where there is a infinite ground plane at the aperture ']);
                    input_info.FF.ap = input('Do you want to consider the aperture reflection? 1, YES, 0, NO: [0] ');
                    
                    if isempty(input_info.FF.ap)
                        input_info.FF.ap = 0;
                    end
                    if (input_info.FF.ap == 1) || (input_info.FF.ap == 0)
                        
                        ap = 0;
                    else
                        disp(['Invalid normalization input, try to input again: ']);
                        ap = 1;
                    end
                end
                
                % Far field domain selection ==============================
                
                
                disp(['Domain of far fields: ']);
                
               
               
  
                    thdomain = 1;
                    phdomain = 1;
                    while thdomain
                        [input_info.FF.th0] = input('Elevation angle (start) in degree: default [-90]: ');
                        if isempty(input_info.FF.th0)
                            input_info.FF.th0 = -90; 
                        end

                        [input_info.FF.thF] = input('Elevation angle (end) in degree: default [90]: ');
                        if isempty(input_info.FF.thF)
                            input_info.FF.thF = 90; 
                        end

                        if input_info.FF.thF > input_info.FF.th0
                            thdomain = 0;
                            Nth = 1;
                            while Nth 
                                [input_info.FF.Nth] = input('Number of elevation points: default [180]: ');
                                if isempty(input_info.FF.Nth)
                                    input_info.FF.Nth = 180;
                                end
                                if (mod(input_info.FF.Nth, 1) == 0) && (input_info.FF.Nth > 0)
                                    Nth = 0;
                                    input_info.FF.th = linspace(input_info.FF.th0, input_info.FF.thF, input_info.FF.Nth);
                                else
                                    disp(['Number of elevation points should be an integer, try to input again: ']);
                                    Nth = 1;
                                end
                            end
                            
                        else
                            disp(['Elevation domian is wrong, makes sure the start angle is less than the end angle: ']);
                            thdomain = 1;
                        end
                    end
                     while phdomain
                        [input_info.FF.ph0] = input('Azimuth angle (start) in degree: default [0]: ');
                        if isempty(input_info.FF.ph0)
                            input_info.FF.ph0 = 0; 
                        end

                        [input_info.FF.phF] = input('Azimuth angle (end) in degree: default [360]: ');
                        if isempty(input_info.FF.phF)
                            input_info.FF.phF = 360; 
                        end

                        if input_info.FF.phF > input_info.FF.ph0
                            phdomain = 0;
                            Nph = 1;
                            while Nph 
                                [input_info.FF.Nph] = input('Number of azimuth points: default [360]: ');
                                if isempty(input_info.FF.Nph)
                                    input_info.FF.Nph = 360;
                                end
                                if (mod(input_info.FF.Nph, 1) == 0) && (input_info.FF.Nph > 0)
                                    Nph = 0;
                                    input_info.FF.ph = linspace(input_info.FF.ph0, input_info.FF.phF, input_info.FF.Nph);
                                else
                                    disp(['Number of azimuth points should be an integer, try to input again: ']);
                                    Nph = 1;
                                end
                            end
                            
                        else
                            disp(['Azimuth domian is wrong, makes sure the start angle is less than the end angle: ']);
                            phdomain = 1;
                        end
                    end
                    
               
                
            else
                disp('Excitation input is wrong, try to input it again: ');
                mff_in = 1; 
                
            end
            
            
       end
 
   
   change_options = input('Give 0 is you agree with this options, 1 if you want to reintroduce them [0]:');
    if isempty(change_options)
        change_options = 0;
    end
    
end


%% =============================== Calculating the Far fields ==============================================


disp('Calling the Far Field function: ');

[Eth, Eph, Eco, Exp, E, Gamma, Dm, R11, Rap] = FF(input_info);

ff.Eth = Eth;
ff.Eph = Eph;
ff.Eco = Eco;
ff.Exp = Exp;
ff.E = E;
ff.f_axis = input_info.FF.f_axis;
ff.th = input_info.FF.th;
ff.ph = input_info.FF.ph;
ff.ex = input_info.FF.ex_mode;
ff.ar = input_info.FF.ar; 
ff.norm = input_info.FF.norm;
ff.mff = input_info.FF.mff;
ff.ap = input_info.FF.ap;
ff.GSM = input_info.FF.GSM;

if input_info.FF.ap == 1
    ff.Gamma = Gamma;
    ff.Dm = Dm;
    ff.R11 = R11;
    ff.Rap = Rap;
end

%% ============================================== Save Data ==============================================================

%% Far Fields

disp("Save Data Section: ")

save_ = input('Do you want to save the results of Far Fields in a file? 1 for Yes, 0 for No: [1]: ');

if isempty(save_)
    save_ = 0;
end

if (save_ == 1)
    disp('Selct the folder in the UI: ')
    selected_dir = uigetdir();
    file_name = input("Input the file name you want: Use single quotes when writing file name such as: 'Test': ");
   
    save([selected_dir, '/', file_name], 'ff');
end

disp(['Data saved in', fullfile(selected_dir,file_name), '.mat']);

disp(['For plotting the Far fiels from a mat file, please use Plot_FF.m File']);

%% end ===================================================================================================================