%% Deliverable#3 Tworit v1.0 S Parameters and Far Fields
%
% THIS SCRIPT IS FOR REFLECTION COEFFICIENTS AND FAR FIELDS OF CONICAL HORN ANTENNA
%
% DESCRIPTION: 
% To find out reflection coefficients and far fields of a cascaded
% cylindrical waveguide structure/ conical horn structure


% SETTINGS:
% Setting ¨a¨ is all about the geometry of the waveguide/ horn and setting
% ¨b¨ is all about the waveguide modes
% ¨c¨ is all about far fields inputs

%        a) GEOMETRY: Input if you want to simulate a cascaded structure of
%        cylindrical waveguides or a cascaded structure of cones 
%        In case of a cascaded structure of cones, the discritization is
%        done based on one tenth of a wavelength (lambda/10) of the maximum frequency
%        provided. For every cylindrical/conical waveguides you need to
%        provide the dimensions of the structures
%        b) Far field configuration [frequency, excitation of modes, azimuth, elevation,
%        normalization] and an option to include/ exclude aperture refletion

% OUTPUT: 
%        a) Reflection/Transmission Coefficients of mode pairs in 3D
%        format:It can be saved with an option

%           DIM1: Frequency axis
%           DIM2: Modes at the base waveguide
%           DIM3: Modes at the top waveguide

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
   ref_ff = 1;
   while ref_ff
    
       input_info.REF_FF = input('Select 1 if you want to compute only the reflection coefficient, 2 if reflection coefficients and far fields [1]:  ');
       if isempty(input_info.REF_FF)
        input_info.REF_FF = 1;
       
       end
       
       if (input_info.REF_FF == 1) || (input_info.REF_FF == 2)
            ref_ff = 0; 
       end
        
       if ref_ff == 1
            disp('This input is invalid, try entering again: ');
       end
       
   end
   
   if (input_info.REF_FF == 1) || (input_info.REF_FF == 2)
%% ==================================================== Geometry of the Antenna ============================================================
    type = 1;
    
    while type 
        input_info.GEO.type = input('Select, 1: Cascaded cylinders, Select, 2: Cylinder base with Cascaded Cone structure: ');
        
        if (input_info.GEO.type == 1) || (input_info.GEO.type == 2)
            type = 0;
        end
        
        if type == 1
            disp('The Geometry input type is invalid, please type either 1 or 2 and enter: ');
        end
        
        if input_info.GEO.type == 1
            N = 1;
            while N
            input_info.GEO.N = input('Number of cylinders you want: ');
            
            if (mod(input_info.GEO.N, 1) == 0) && ( input_info.GEO.N > 0 )
                N = 0;
            else
                disp('Invalid number of cylinders, try to input again: ');
                N = 1;
            end
            end
            
            for i = 1:input_info.GEO.N
                input_info.GEO.E(i).r = input(['Radius of ', ' Element ', num2str(i) , ' in [m]: ']);
                input_info.GEO.E(i).l = input(['Length of ', ' Element ', num2str(i) , ' in [m]: ']); 
                input_info.GEO.E(i).er = input(['Relative permittivity of  ', ' Element ', num2str(i) , ' [1]: ']);
                 if isempty(input_info.GEO.E(i).er )
                    input_info.GEO.E(i).er = 1;
                 end
                input_info.GEO.E(i).mur = input(['Relative permeability of  ', ' Element ', num2str(i) , ' [1]: ']);
                if isempty(input_info.GEO.E(i).mur )
                    input_info.GEO.E(i).mur = 1;
                end
            end
        elseif input_info.GEO.type == 2
               N = 1;
               while N
                    input_info.GEO.N = input('Number of cones you want: ');
            
                    if (mod(input_info.GEO.N, 1) == 0) && ( input_info.GEO.N > 0 )
                         N = 0;
                    else
                        disp('Invalid number of cones, try to input again: ');
                        N = 1;
                    end
                end
            

                input_info.GEO.E(1).r = input(['Radius of the base cylinder Element ', ' in [m]: ']);
                input_info.GEO.E(1).l = input(['Length of the base cylinder Element ', 'in [m]: ']);
                input_info.GEO.E(1).er = input(['Relative permittivity of the base cylinderin [1]: ']);
                if isempty(input_info.GEO.E(1).er )
                    input_info.GEO.E(1).er = 1;
                 end
                input_info.GEO.E(1).mur = input(['Relative permeability of the base cylinderin [1]: ']);
                if isempty(input_info.GEO.E(1).mur )
                    input_info.GEO.E(1).mur = 1;
                end
                
            for i = 1:input_info.GEO.N
                disp('The base radius of every conical element is same as the radius of the previous geometrical element: ');
                disp('========================================================================================');
                input_info.GEO.E(i+1).r = input(['Radius of the top of ', 'Conical Element ', num2str(i) ,  ' in [m]: ']);
                input_info.GEO.E(i+1).l = input(['Length of ', 'Conical Element ', num2str(i) ,  ' in [m]: ']);
                input_info.GEO.E(i+1).er = input(['Relative Permittivity of ', 'Conical Element ', num2str(i) ,  ' [1]: ']);
                if isempty(input_info.GEO.E(i+1).er )
                    input_info.GEO.E(i+1).er = 1;
                 end
                input_info.GEO.E(i+1).mur = input(['Relative Permittivity of ', 'Conical Element ', num2str(i) ,  ' [1]: ']);
                if isempty(input_info.GEO.E(i+1).mur )
                    input_info.GEO.E(i+1).mur = 1;
                end
            end
        end
        
        
    end
    
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

        
   end
   
   if (input_info.REF_FF == 2)
        
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
                disp(['The excitation chosen is the mode ']);
                for mffi = 1:length(mff)
                    disp(modes(mffi).mode); disp(modes(mffi).m); disp(modes(mffi).n); disp('Orientation: '); disp(modes(mffi).pol * 180/pi); disp(' [deg]'); 
                end
                input_info.FF.ex_mode = modes(mff);
                input_info.FF.mff = mff;
                
                % Frequency Selection =====================================
                disp(['Frequencies should be in between ', num2str(f_cutoff(min(mff))*1e-9), 'GHz and ', num2str(input_info.MODES.f_axis(end)*1e-9), 'GHz']);
                
                fopt = 1;
                while fopt 
                
                    input_info.FF.fopt = input('Do you want to consider all the frequency points in this range? 1 YES, 0 No [1]:');
                    if isempty(input_info.FF.fopt)
                        input_info.FF.fopt = 1;
                    end
                    
                    if (input_info.FF.fopt == 0) || (input_info.FF.fopt == 1)
                        fopt = 0;
                        [~, fidx] = min(abs(input_info.MODES.f_axis - f_cutoff(min(mff))));
                        if input_info.MODES.f_axis(fidx) < f_cutoff(min(mff))
                            fidx = fidx + 1;
                        end
                        
                        f_axis_ff = input_info.MODES.f_axis(fidx:end);
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

%% =============================== Calculating the Far fields ==============================================
if input_info.REF_FF == 2

disp('Calling the Far Field function: ');


ff.SRR = b.SRR(input_info.FF.Nfi, :, :);
ff.SRT = b.SRT(input_info.FF.Nfi, :, :);
ff.STR = b.STR(input_info.FF.Nfi, :, :);
ff.STT = b.STT(input_info.FF.Nfi, :, :);

input_info.FF.GSM = ff;

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

end



%% ============================================== Save Data ==============================================================
%% Geometrical data
geo = input_info.GEO;


disp("Save Data Section: ")

save_ = input('Do you want to save the Geometry in a file? 1 for Yes, 0 for No: [1]: ');

if isempty(save_)
    save_ = 0;
end

if (save_ == 1)
    disp('Selct the folder in the UI: ')
    selected_dir = uigetdir();
    file_name = input("Input the file name you want: Use single quotes when writing file name such as: 'Test': ");
   
    save([selected_dir, '/', file_name], 'geo');
    disp(['Data saved in', fullfile(selected_dir,file_name), '.mat']);
end



disp(['For plotting the Geometry from a mat file, please use Plot_Geomat.m File']);

%% S parameters

disp("Save Data Section: ")

save_ = input('Do you want to save the results of S parameters in a file? 1 for Yes, 0 for No: [1]: ');

if isempty(save_)
    save_ = 0;
end

if (save_ == 1)
    disp('Selct the folder in the UI: ')
    selected_dir = uigetdir();
    file_name = input("Input the file name you want: Use single quotes when writing file name such as: 'Test': ");
   
    save([selected_dir, '/', file_name], 'b');
    
    disp(['Data saved in', fullfile(selected_dir,file_name), '.mat']);
end


disp(['For plotting the S parameters from a mat file, please use Plot_Sparam.m File']);

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
    disp(['Data saved in', fullfile(selected_dir,file_name), '.mat']);
end



disp(['For plotting the Far fields from a mat file, please use Plot_FF.m File']);

%% end ===================================================================================================================

