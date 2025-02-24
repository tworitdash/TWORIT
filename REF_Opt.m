%% Deliverable#9 Tworit v3.0 S Parameter Optimization
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

%        a) GEOMETRY: Input if you want to simulate a cascaded structure of
%        cylindrical waveguides or a cascaded structure of cones 
%        In case of a cascaded structure of cones, the discritization is
%        done based on one tenth of a wavelength (lambda/10) of the maximum frequency
%        provided. For every cylindrical/conical waveguides you need to
%        provide the dimensions of the structures.
%           

% OUTPUT: 
%        a) Optimized dimensions of the antenna - Complete Geometry

%% 
clear;
close all;
addpath('./Tworit_functions_toolbox/');
change_options = 1;
while change_options
    
    
%% ==================================================== Geometry of the Antenna ============================================================
    type = 1;
    
    while type 
        
        disp('===================== Input the design space for the Geometry ===================================== ');
        input_info.GEO.type = input('Select, 1: Cascaded cylinders, Select, 2: Cylinder base with Cascaded Cone structure: ');
        
        if (input_info.GEO.type == 1) || (input_info.GEO.type == 2)
            type = 0;
        end
        
        if type == 1
            disp('The Geometry input type is invalid, please type either 1 or 2 and enter: ');
        end
        
        if input_info.GEO.type == 1
            input_info.GEO.N = input('Number of cylinders you want: ');
            for i = 1:input_info.GEO.N
                input_info.GEO.E(i).r_l = input(['Lower Limit for the Radius of ', ' Element ', num2str(i) , ' in [m]: ']);
                input_info.GEO.E(i).r_h = input(['Higher Limit for the Radius of ', ' Element ', num2str(i) , ' in [m]: ']);
                input_info.GEO.E(i).r_0 = input(['Start point for optimization for the Radius of ', ' Element ', num2str(i) , ' in [m]: ']);
                
                input_info.GEO.E(i).l_l = input(['Lower Limit for the Length of ', ' Element ', num2str(i) , ' in [m]: ']);
                input_info.GEO.E(i).l_h = input(['Higher Limit for the Length of ', ' Element ', num2str(i) , ' in [m]: ']);
                input_info.GEO.E(i).l_0 = input(['Start point for optimization for the Length of ', ' Element ', num2str(i) , ' in [m]: ']);
                
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
            input_info.GEO.N = input('Number of cones you want: ');

                input_info.GEO.E(1).r_l = input(['Lower Limit of the Radius of the base cylinder Element ', ' in [m]: ']);
                input_info.GEO.E(1).r_h = input(['Higher Limit of the Radius of the base cylinder Element ', ' in [m]: ']);
                input_info.GEO.E(1).r_0 = input(['Start point for optimization of the Radius of the base cylinder Element ', ' in [m]: ']);
                
                
                input_info.GEO.E(1).l_l = input(['Lower Limit of Length of the base cylinder Element ', 'in [m]: ']);
                input_info.GEO.E(1).l_h = input(['Higher Limit of Length of the base cylinder Element ', 'in [m]: ']);
                input_info.GEO.E(1).l_0 = input(['Start point for the optimization of Length of the base cylinder Element ', 'in [m]: ']);
                
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
                input_info.GEO.E(i+1).r_l = input(['Lower Limit of the Radius of the top of ', 'Conical Element ', num2str(i) ,  ' in [m]: ']);
                input_info.GEO.E(i+1).r_h = input(['Higher Limit of the Radius of the top of ', 'Conical Element ', num2str(i) ,  ' in [m]: ']);
                input_info.GEO.E(i+1).r_0 = input(['Start point for the optimization of the Radius of the top of ', 'Conical Element ', num2str(i) ,  ' in [m]: ']);
                
                input_info.GEO.E(i+1).l_l = input(['Lower Limit of the Length of ', 'Conical Element ', num2str(i) ,  ' in [m]: ']);
                input_info.GEO.E(i+1).l_h = input(['Higher Limit of the Length of ', 'Conical Element ', num2str(i) ,  ' in [m]: ']);
                input_info.GEO.E(i+1).l_0 = input(['Start point for the optimization of the Length of ', 'Conical Element ', num2str(i) ,  ' in [m]: ']);
                
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
%% Just for plotting of the Geometry

    input_info.Plot_Geometry = input('Do you want to plot the Geometry?: Select 1 for YES, Select 2 for NO, Default is NO: ');
    
    if isempty(input_info.Plot_Geometry)
        input_info.Plot_Geometry = 0;
    end
    
    if input_info.Plot_Geometry == 1
        input_info.Plot_Geometry_sub = input('Which configuration to plot?: Select 1 for Lower Limits, Select 2 for Higher Limits, Default is 1: ');
        if isempty(input_info.Plot_Geometry_sub) || (input_info.Plot_Geometry_sub ~= 1) || (input_info.Plot_Geometry_sub ~= 2)
            input_info.Plot_Geometry_sub = 1;
        end
        if input_info.Plot_Geometry_sub == 1
            if input_info.GEO.type == 1
              for o = 1:input_info.GEO.N
                 input_info.GEO.E(o).r = input_info.GEO.E(o).r_l;
                 input_info.GEO.E(o).l = input_info.GEO.E(o).l_l;
              end
              Plot_GEO(input_info.GEO);
            else
              for o = 1:input_info.GEO.N+1
                 input_info.GEO.E(o).r = input_info.GEO.E(o).r_l;
                 input_info.GEO.E(o).l = input_info.GEO.E(o).l_l;
              end
              Plot_GEO(input_info.GEO);
            end
        else
            if input_info.GEO.type == 1
              for o = 1:input_info.GEO.N
                 input_info.GEO.E(o).r = input_info.GEO.E(o).r_h;
                 input_info.GEO.E(o).l = input_info.GEO.E(o).l_h;
              end
              Plot_GEO(input_info.GEO);
            else
              for o = 1:input_info.GEO.N+1
                 input_info.GEO.E(o).r = input_info.GEO.E(o).r_h;
                 input_info.GEO.E(o).l = input_info.GEO.E(o).l_h;
              end
              Plot_GEO(input_info.GEO);
            end
            
        end
    end
    
%% ========================================== Frequency and modes ==============================================

f_invalid = 1;
    
%     while f_invalid
    input_info.MODES.f_begin = input('Input the frequency of the optimization in [GHz]: ')*1e9;

    change_options = input('Give 0 if you agree with this options, 1 if you want to reintroduce them [0]:');
    if isempty(change_options)
        change_options = 0;
    end
end

%% =========================== Optimizer code ==================================================================

input_info.Opt.type = input('Which Optimizer you want to use: 1 for FMINCON, 2 for GENETIC ALGORITHM: [1]');

if isempty(input_info.Opt.type) % || (input_info.Opt.type ~= 1) || ( input_info.Opt.type ~= 2 )
    input_info.Opt.type = 1;
end


if input_info.Opt.type == 1
    
    Opt.solver = 'fmincon';
    
    if input_info.GEO.type == 1
        for k = 1:input_info.GEO.N
            lb(k) = input_info.GEO.E(k).r_l;
            ub(k) = input_info.GEO.E(k).r_h;
            x0(k) = input_info.GEO.E(k).r_0;
        end
        for k = input_info.GEO.N+1:input_info.GEO.N*2
            lb(k) = input_info.GEO.E(k-input_info.GEO.N).l_l;
            ub(k) = input_info.GEO.E(k-input_info.GEO.N).l_h;
            x0(k) = input_info.GEO.E(k_input_info.GEO.N).l_0;
        end
        
        Opt.objective = @(x) GSM_v1_OPT(x(1:2*input_info.GEO.N), input_info);
    else
        for k = 1:input_info.GEO.N+1
            lb(k) = input_info.GEO.E(k).r_l;
            ub(k) = input_info.GEO.E(k).r_h;
            x0(k) = input_info.GEO.E(k).r_0;
        end
        for k = input_info.GEO.N+2:(input_info.GEO.N+1)*2
            lb(k) = input_info.GEO.E(k-input_info.GEO.N-1).l_l;
            ub(k) = input_info.GEO.E(k-input_info.GEO.N-1).l_h;
            x0(k) = input_info.GEO.E(k-input_info.GEO.N-1).l_0;
        end
        Opt.objective = @(x) GSM_v1_OPT(x(1:2*(input_info.GEO.N+1)), input_info);
    end
    
    Opt.lb = [lb.'];
    Opt.ub = [ub.'];
    Opt.x0 = [x0.'];
    
    Opt.options = optimoptions(@fmincon, 'PlotFcn', ...
        {'optimplotx', 'optimplotfirstorderopt'}, 'UseParallel', true, 'MaxIterations', 20);
    tic;
    [Result, fval, exf2, out] = fmincon(Opt);
    time_consumed = toc;
    
else
    
    Opt.solver = 'ga';
    
    if input_info.GEO.type == 1
        for k = 1:input_info.GEO.N
            lb(k) = input_info.GEO.E(k).r_l;
            ub(k) = input_info.GEO.E(k).r_h;
            x0(k) = input_info.GEO.E(k).r_0;
            
        end
        for k = input_info.GEO.N+1:input_info.GEO.N*2
            lb(k) = input_info.GEO.E(k-input_info.GEO.N).l_l;
            ub(k) = input_info.GEO.E(k-input_info.GEO.N).l_h;
            x0(k) = input_info.GEO.E(k_input_info.GEO.N).l_0;
            
        end
        
        nvars = 2*input_info.GEO.N;
        
        Opt.fitnessfcn = @(x) GSM_v1_OPT(x(1:2*input_info.GEO.N), input_info);
    else
        for k = 1:input_info.GEO.N+1
            lb(k) = input_info.GEO.E(k).r_l;
            ub(k) = input_info.GEO.E(k).r_h;
            x0(k) = input_info.GEO.E(k).r_0;
        end
        for k = input_info.GEO.N+2:(input_info.GEO.N+1)*2
            lb(k) = input_info.GEO.E(k-input_info.GEO.N-1).l_l;
            ub(k) = input_info.GEO.E(k-input_info.GEO.N-1).l_h;
            x0(k) = input_info.GEO.E(k-input_info.GEO.N-1).l_0;
        end
        Opt.nvars = 2*(input_info.GEO.N+1);
        Opt.fitnessfcn = @(x) GSM_v1_OPT(x(1:2*(input_info.GEO.N+1)), input_info);
    end
    
    Opt.lb = [lb.'];
    Opt.ub = [ub.'];
    IP = [x0];
    
    Opt.fitnesslimit = input('Input the minimum limit till which the solver needs to reach in dB: ');
    
    Opt.options = optimoptions(@ga, 'PlotFcn', {'gaplotbestf', 'gaplotbestindiv'}, 'Display', 'iter',... 
    'InitialPopulationMatrix', [IP], 'UseParallel',...
    true, 'FitnessLimit', Opt.fitnesslimit );
    tic;
    [Result, fval, exf2, out] = ga(Opt);
    time_consumed = toc;
end

%% Save Data Section

output_info.GEO.type = input_info.GEO.type;
output_info.GEO.N = input_info.GEO.N;

if output_info.GEO.type == 1
    for p = 1:output_info.GEO.N
        output_info.GEO.E(p).r = Result(p);
        output_info.GEO.E(p).l = Result(output_info.GEO.N+p);
        output_info.GEO.E(p).er = input_info.GEO.E(p).er;
        output_info.GEO.E(p).mur = input_info.GEO.E(p).mur;
    end
else
    for p = 1:output_info.GEO.N+1
        output_info.GEO.E(p).r = Result(p);
        output_info.GEO.E(p).l = Result(output_info.GEO.N+1+p);
        output_info.GEO.E(p).er = input_info.GEO.E(p).er;
        output_info.GEO.E(p).mur = input_info.GEO.E(p).mur;
    end
end

geo = output_info.GEO;

disp("Save Data Section: ")

save_ = input('Do you want to save the geometrical results in a file? 1 for Yes, 0 for No: [1]: ');

if isempty(save_)
    save_ = 0;
end

if (save_ == 1)
    disp('Selct the folder in the UI: ')
    selected_dir = uigetdir();
    file_name = input("Input the file name you want: Use single quotes when writing file name such as: 'Test': ");
   
    save([selected_dir, '/', file_name], 'geo');
end

%% Save Optimiztion related outputs

optim.results = Result;
optim.fval = fval;
optim.exf2 = exf2;
optim.out = out;
optim.time_consumed = time_consumed;

disp("Save Data Section: ")

save_ = input('Do you want to save the optimization result configuration in a file? 1 for Yes, 0 for No: [1]: ');

if isempty(save_)
    save_ = 0;
end

if (save_ == 1)
    disp('Selct the folder in the UI: ')
    selected_dir = uigetdir();
    file_name = input("Input the file name you want: Use single quotes when writing file name such as: 'Test': ");
   
    save([selected_dir, '/', file_name], 'optim');
end

%% Save optimization related inputs

optim.results = Result;
optim.fval = fval;
optim.exf2 = exf2;
optim.out = out;

disp("Save Data Section: ")

save_ = input('Do you want to save the optimization input configuration in a file? 1 for Yes, 0 for No: [1]: ');

if isempty(save_)
    save_ = 0;
end

if (save_ == 1)
    disp('Selct the folder in the UI: ')
    selected_dir = uigetdir();
    file_name = input("Input the file name you want: Use single quotes when writing file name such as: 'Test': ");
   
    save([selected_dir, '/', file_name], 'Opt');
end


%% Save Data section


% disp('Calling the S paremeter function: ');
% 
% 
% 
% 
% % [RLRR_TE11, Srr, Stt, Str, Srt] = GSM_v1(input_info);
% 
% 
% 
% 
% disp('S Parameters Computed: ');
% 
% 
% % figure; plot(input_info.MODES.f_axis .* 1e-9, db(STT(:, 1, 1)));
% 
% b.STT = Stt;
% b.SRR = Srr;
% b.STR = Str;
% b.SRT = Srt;
% b.RLRR_TE11 = RLRR_TE11;
% b.f_axis = input_info.MODES.f_axis;
% 
% %% ============================================== Save Data ==============================================================
% %% Geometrical data
% geo = input_info.GEO;
% 
% 
% disp("Save Data Section: ")
% 
% save_ = input('Do you want to save the Geometry in a file? 1 for Yes, 0 for No: [1]: ');
% 
% if isempty(save_)
%     save_ = 0;
% end
% 
% if (save_ == 1)
%     disp('Selct the folder in the UI: ')
%     selected_dir = uigetdir();
%     file_name = input("Input the file name you want: Use single quotes when writing file name such as: 'Test': ");
%    
%     save([selected_dir, '/', file_name], 'geo');
% end
% 
% disp(['Data saved in', fullfile(selected_dir,file_name), '.mat']);
% 
% disp(['For plotting the Geometry from a mat file, please use Plot_Geomat.m File']);
% 
% %% S Parameters
% 
% disp("Save Data Section: ")
% 
% save_ = input('Do you want to save the results in a file? 1 for Yes, 0 for No: [1]: ');
% 
% if isempty(save_)
%     save_ = 0;
% end
% 
% 
% 
% if (save_ == 1)
%     disp('Selct the folder in the UI: ')
%     selected_dir = uigetdir();
%     file_name = input("Input the file name you want: Use single quotes when writing file name such as: 'Test': ");
%    
%     save([selected_dir, '/', file_name], 'b');
% end
% 
% 
% disp(['Data saved in', fullfile(selected_dir,file_name), '.mat']);
