addpath('./Tworit_functions_toolbox/');
%% Deliverable#5 Tworit v1.0 S Plotting Far Fields from saved Far Fields parameters in .mat files

% Inputs appear as the code is run:

%% ==================================================== Plot Section =======================================================

Plot_enable = input('Do you want to plot the Far Fields you like now? : 1 for Yes, 0 for No: [0]: ');

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

FF = load(plot_info.file);

FF = FF.ff;

% disp("T is the top waveguide: R is the base Waveguide: ")
% 
disp("Select the independent variable: ")

FF_valid = 1;

while FF_valid

  Config_FF = input('1 for frequency, 2 for elevation, 3 for azimuth: Default is [2]: ');
  
  if isempty(Config_FF)
      Config_FF = 2;
      
  end
  
  if (Config_FF == 1) || (Config_FF == 2) || (Config_FF == 3)
      FF_valid  = 0;
      if Config_FF == 1
          
          idf = 1:length(FF.f_axis);
          x = FF.f_axis .* 1e-9;
          ele = 1;
          while ele
              disp('Choose an index for the elevation angle axis: The elevation angles and their indices are shown below');
              for i = 1:length(FF.th)
                    disp('index: '); disp(i);
                    disp('index: '); disp(FF.th(i));
              end
              idth = input('Choose an index from the above output for elevation: default [1]: ');
              if isempty(idth)
                  idth = 1;
              end
              if (idth < 1) || (idth > length(FF.th))
                  disp(['Wrong input for elevation index: try one index from 1 to ', num2str(length(FF.th))]);
                  ele = 1;
              else
                  ele = 0;
              end
          end
          
          azi = 1;
          while azi
              disp('Choose an index for the azimuth angle axis: The azimuth angles and their indices are shown below');
              for i = 1:length(FF.ph)
                    disp('index: '); disp(i);
                    disp('index: '); disp(FF.ph(i));
              end
              idph = input('Choose an index from the above output for azimuth: default [1]: ');
              if isempty(idph)
                  idph = 1;
              end
              if (idph < 1) || (idph > length(FF.ph))
                  disp(['Wrong input for azimuth index: try one index from 1 to ', num2str(length(FF.ph))]);
                  azi = 1;
              else
                  azi = 0;
              end
          end
          
          lgd = ['FF at \theta  = ', num2str(FF.th(idth)), ' [\circ] and \phi  = ', num2str(FF.ph(idph)), ' [\circ]'];
          xl = ['Frequency [GHz] '];
       elseif Config_FF == 2
           
          idth = 1:length(FF.th);
          x = FF.th;
          azi = 1;
          while azi
              disp('Choose an index for the azimuth angle axis: The azimuth angles and their indices are shown below');
              for i = 1:length(FF.ph)
                    disp('index: '); disp(i);
                    disp('index: '); disp(FF.ph(i));
              end
              idph = input('Choose an index from the above output for azimuth: default [1]: ');
              if isempty(idph)
                  idph = 1;
              end
              if (idph < 1) || (idph > length(FF.ph))
                  disp(['Wrong input for azimuth index: try one index from 1 to ', num2str(length(FF.ph))]);
                  azi = 1;
              else
                  azi = 0;
              end
          end
          
          fre = 1;
          while fre
              disp('Choose an index for the frequency axis: The frequencies and their indices are shown below');
              for i = 1:length(FF.f_axis)
                    disp('index: '); disp(i);
                    disp('index: '); disp(FF.f_axis(i)*1e-9); disp(' GHz');
              end
              idf = input('Choose an index from the above output for frequencies: default [1]: ');
              if isempty(idf)
                  idf = 1;
              end
              if (idf < 1) || (idf > length(FF.f_axis))
                  disp(['Wrong input for frequency index: try one index from 1 to ', num2str(length(FF.f_axis))]);
                  fre = 1;
              else
                  fre = 0;
              end
          end
          
          lgd = ['FF at F  = ', num2str(FF.f_axis(idf).*1e-9), ' [GHz] and \phi  = ', num2str(FF.ph(idph)), ' [\circ]'];
          xl = ['Elevation [\circ] '];
           
      else
          idph = 1:length(FF.ph);
          x = FF.ph;
          ele = 1;
          while ele
              disp('Choose an index for the elevation angle axis: The elevation angles and their indices are shown below');
              for i = 1:length(FF.th)
                    disp('index: '); disp(i);
                    disp('index: '); disp(FF.th(i));
              end
              idth = input('Choose an index from the above output for elevation: default [1]: ');
              if isempty(idth)
                  idth = 1;
              end
              if (idth < 1) || (idth > length(FF.th))
                  disp(['Wrong input for elevation index: try one index from 1 to ', num2str(length(FF.th))]);
                  ele = 1;
              else
                  ele = 0;
              end
          end
          
          fre = 1;
          while fre
              disp('Choose an index for the frequency axis: The frequencies and their indices are shown below');
              for i = 1:length(FF.f_axis)
                    disp('index: '); disp(i);
                    disp('index: '); disp(FF.f_axis(i)*1e-9); disp(' GHz');
              end
              idf = input('Choose an index from the above output for frequencies: default [1]: ');
              if isempty(idf)
                  idf = 1;
              end
              if (idf < 1) || (idf > length(FF.f_axis))
                  disp(['Wrong input for frequency index: try one index from 1 to ', num2str(length(FF.f_axis))]);
                  fre = 1;
              else
                  fre = 0;
              end
          end
          lgd = ['FF at F  = ', num2str(FF.f_axis(idf).*1e-9), ' [GHz] and \phi  = ', num2str(FF.ph(idph)), ' [\circ]'];
          xl = ['Elevation [\circ] '];
       end
   end
  
  if FF_valid == 1
      disp('Given configuration is wrong: Please input one of the following values only: ');
  end
  
end

ff_conf = 1;

while ff_conf

out = input('What do you want to plot?, 1 for Eth, 2 for Eph, 3 for E, 4 for Eco, 5 for Exp: Default is [3]: ');

if isempty(out)
    out = 3;
end

if (out == 1) || (out == 2) || (out == 3) || (out == 4) || (out == 5)
    ff_conf = 0;
    if out == 1
        y = db(squeeze(FF.Eth(idf, idph, idth)));
        yl = ['E Field [dB]'];
        lgd = ['E_{\theta} [dB]', lgd];
    elseif out == 2
        y = db(squeeze(FF.Eph(idf, idph, idth))); 
        lgd = ['E Field [dB]'];
        
        yl = ['E_{\phi} [dB]', lgd];
    elseif out == 3
        y = db(squeeze(FF.E(idf, idph, idth)));
        yl = ['E Field [dB]'];
        lgd = ['E_{total} [dB]', lgd];
    elseif out == 4
        y = db(squeeze(FF.Eco(idf, idph, idth)));
        yl = ['E Field [dB]'];
        lgd = ['E_{CO} [dB]', lgd];
    else
        y = db(squeeze(FF.Exp(idf, idph, idth)));
        yl = ['E Field [dB]'];
        lgd = ['E_{XP} [dB]', lgd];
    end
else
    ff_conf = 1;
    disp("Wrong configuration, please input again: ")
end
end

plot_info.x = x;
plot_info.y = y;
plot_info.xl = xl;
plot_info.yl = yl;
plot_info.lgd = lgd;


plot_FF(plot_info);

Plot_enable = input('Want to plot a new Far Field?: 1 for Yes and 0 for No: [0]: ');
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