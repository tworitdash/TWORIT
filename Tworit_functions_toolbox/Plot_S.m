function [] = Plot_S(plot_info)

%% Plot_S plots the S parameters from a saved file with S parameter configuration of a cascaded cylindrical waveguide/ conical waveguide problems.

% Inputs: 
%       plot_info:            It is a structure input having the following fileds
%       plot_info.file -      Full file path to the .mat file containing the S
%                             parameters. Check REF.m to see how this file is
%                             created fromt this solver.
%       plot_info.Config_S  - S parameter configuration: 1 for SRR, 2 for
%                             STT, 3 for STR, 4 for SRT; T is the bigger
%                             and R is the smaller waveguide ends.
%       plot_info.Config_m -  Serial number of the first mode in the
%                             desired S parameter plot. To find the serial
%                             number of a desired mode, check MODES.csv in
%                             the parent directory of the tool
%       plot_info.Config_n -  Serial number of the second mode in the
%                             desired S parameter plot. To find the serial
%                             number of a desired mode, check MODES.csv in
%                             the parent directory of the tool
%       plot_info.same_diff - Plot in the same plot or do you want a
%                             different plot. 1 for same plot, 0 for
%                             different plot
% Outputs: 
%        Doesn't return anything. 
%        Only plots the S parameters in dB for the waveguide pair/ mode pair selected by the user.
% For a proper usage, please check Plot_Sparam.m in the parent directory of the
% tool
%% 
  Str = load('Xmn.mat');
  str = Str.Xmn;
  MM = load(plot_info.file);
  plot_info.f_axis = MM.b.f_axis;
  
  if plot_info.Config_S == 1
     S = MM.b.SRR;
     txt = ['SRR of ', str(plot_info.Config_m).mode, num2str(str(plot_info.Config_m).m), ...
         num2str(str(plot_info.Config_m).n), ' Pol: ', num2str(str(plot_info.Config_m).pol * 180/pi),...
         ' and ', str(plot_info.Config_n).mode, ...
         num2str(str(plot_info.Config_n).m), ...
         num2str(str(plot_info.Config_n).n), ' Pol: ', num2str(str(plot_info.Config_n).pol .* 180/pi)]; 
     txt = char(txt);
     
     if plot_info.same_diff == 0
        figure; plot(plot_info.f_axis .* 1e-9, db(abs(S(:, plot_info.Config_m, plot_info.Config_n))), 'LineWidth', 2, 'DisplayName', txt);
     else
        hold on; plot(plot_info.f_axis .* 1e-9, db(abs(S(:, plot_info.Config_m, plot_info.Config_n))) , 'LineWidth', 2, 'DisplayName', txt);
     end
     
     grid on;
     legend;
     
  elseif plot_info.Config_S == 2 
      
     S = MM.b.STT;
     txt = ['STT of ', str(plot_info.Config_m).mode, num2str(str(plot_info.Config_m).m), ...
         num2str(str(plot_info.Config_m).n), ' Pol: ', num2str(str(plot_info.Config_m).pol * 180/pi),...
         ' and ', str(plot_info.Config_n).mode, ...
         num2str(str(plot_info.Config_n).m), ...
         num2str(str(plot_info.Config_n).n), ' Pol: ', num2str(str(plot_info.Config_n).pol .* 180/pi)]; 
     txt = char(txt);
     if plot_info.same_diff == 0
        figure; plot(plot_info.f_axis .* 1e-9, db(abs(S(:, plot_info.Config_m, plot_info.Config_n))), 'LineWidth', 2, 'DisplayName', txt);
     else
        hold on; plot(plot_info.f_axis .* 1e-9, db(abs(S(:, plot_info.Config_m, plot_info.Config_n))) , 'LineWidth', 2, 'DisplayName', txt);
     end
     
     grid on;
     legend;
      
  elseif plot_info.Config_S == 3
      
      S = MM.b.STR;
    txt = ['STR of ', str(plot_info.Config_m).mode, num2str(str(plot_info.Config_m).m), ...
         num2str(str(plot_info.Config_m).n), ' Pol: ', num2str(str(plot_info.Config_m).pol * 180/pi),...
         ' and ', str(plot_info.Config_n).mode, ...
         num2str(str(plot_info.Config_n).m), ...
         num2str(str(plot_info.Config_n).n), ' Pol: ', num2str(str(plot_info.Config_n).pol .* 180/pi)]; 
     txt = char(txt);
     
     if plot_info.same_diff == 0
        figure; plot(plot_info.f_axis .* 1e-9, db(abs(S(:, plot_info.Config_m, plot_info.Config_n))), 'LineWidth', 2, 'DisplayName', txt);
     else
        hold on; plot(plot_info.f_axis .* 1e-9, db(abs(S(:, plot_info.Config_m, plot_info.Config_n))) , 'LineWidth', 2, 'DisplayName', txt);
     end
     
     grid on;
     legend;
  else 
      
      S = MM.b.SRT;
     txt = ['SRT of ', str(plot_info.Config_m).mode, num2str(str(plot_info.Config_m).m), ...
         num2str(str(plot_info.Config_m).n), ' Pol: ', num2str(str(plot_info.Config_m).pol * 180/pi),...
         ' and ', str(plot_info.Config_n).mode, ...
         num2str(str(plot_info.Config_n).m), ...
         num2str(str(plot_info.Config_n).n), ' Pol: ', num2str(str(plot_info.Config_n).pol .* 180/pi)]; 
     txt = char(txt);
     
     if plot_info.same_diff == 0
        figure; plot(plot_info.f_axis .* 1e-9, db(abs(S(:, plot_info.Config_m, plot_info.Config_n))), 'LineWidth', 2, 'DisplayName', txt);
     else
        hold on; plot(plot_info.f_axis .* 1e-9, db(abs(S(:, plot_info.Config_m, plot_info.Config_n))) , 'LineWidth', 2, 'DisplayName', txt);
     end
     
     grid on;
     legend;
      
   
  end
  
grid on;
xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('S parameters', 'FontSize', 12, 'FontWeight', 'bold');


end