%% Deliverable#7 Tworit v1.0 Form the Geometry of the Antenna
%
% THIS SCRIPT IS FOR CONSTRUCTING THE GEOMETRY OF THE ANTENNA
%
% DESCRIPTION: 
% To find out reflection coefficients of a cascaded
% cylindrical waveguide structure/ conical horn structure


% SETTINGS:

%        a) GEOMETRY: Input if you want to simulate a cascaded structure of
%        cylindrical waveguides or a cascaded structure of cones 
%        In case of a cascaded structure of cones, the discritization is
%        done based on one tenth of a wavelength (lambda/10) of the maximum frequency
%        provided. For every cylindrical/conical waveguides you need to
%        provide the dimensions of the structures.

%       
%           

% OUTPUT: 
%        a) Geometry with all the radii and length of the elements
%% 
clear;
close all;

addpath('./Tworit_functions_toolbox/');
change_options = 1;
while change_options
    
    
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
            input_info.GEO.N = input('Number of cylinders you want: ');
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
            input_info.GEO.N = input('Number of cones you want: ');

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
    
    change_options = input('Give 0 is you agree with this options, 1 if you want to reintroduce them [0]:');
    if isempty(change_options)
        change_options = 0;
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
end

disp(['Data saved in', fullfile(selected_dir,file_name), '.mat']);

disp(['For plotting the Geometry from a mat file, please use Plot_Geomat.m File']);



%% end ===================================================================================================================

