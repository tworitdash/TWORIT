function [RLRR_TE11, SRR, STT, STR, SRT] = GSM_v1(input_info)

%% GSM_v1 returns the S parameters and the return loss based on the fundamental 
% mode excitaion for a cascaded cylindrical waveguide/ conical waveguide with a 
% cylindrical waveguide at the base
% Inputs: 
% Geometrical Inputs:
%       input_info:               It is a structure input having the following fileds
%       input_info.GEO.type       - 1 for cascaded cylindrical waveguide structure, 2 for
%                                 cascaded conical waveguide structure with a cylindrical waveguide
%                                 as the base
%       input_info.GEO.N          - If GEO.type is 1, it is the number of the cylindrical
%                                 elements. 
%                                 If GEO.type is 2, it is the number of conical
%                                 elements. So for a conical waveguide structure, the
%                                 number of array elements in GEO.E() is N + 1 as it also
%                                 has to contain the dimensions of the base waveguide.
%       input_info.GEO.E()        - Structure Array with the following fields:
%       input_info.GEO.E().r      - Radius of every element
%       input_info.GEO.E().l      - Length of every element
%       input_info.GEO.E(i).er    - Relative Permittivity of the waveguide
%       input_info.GEO.E(i).mur   - Relative Permeability of the waveguide
% Modal Inputs:
%       input_info.MODES.f_axis - vector with frequencies of operation
%       input_info.MODES.Nbn    - Number of modes in the base waveguide
%       input_info.MODES.Ntn    - Number of modes in the top waveguide
% Outputs:
%       RLRR_TE11               - Return loss in dB with a fundamental mode
%                                 excitation at the base waveguide
%       [SRR, STT, STR, SRT]    - GSM of the entire sturcture - R is the
%                                 base waveguide and T is the top waveguide
% Check REF.m in the parent foler of the tool for full implementation

disp('====================================================================================================');


c0 = 3e8;

if input_info.GEO.type == 1
    for i = 1:input_info.GEO.N
        R(i) = input_info.GEO.E(i).r;
        L(i) = input_info.GEO.E(i).l;
        er(i) = input_info.GEO.E(i).er;
        mur(i) = input_info.GEO.E(i).mur;
    end 
else
   disp('Calculating the Number of discrete geometrical elements in the Conical structure based on 1/10 of the minimum wavelength chosen: ');
%    lamb_opt_freq = c0./input_info.MODES.f_axis;
   
   R(1) = input_info.GEO.E(1).r;
   L(1) = input_info.GEO.E(1).l;
   
   er(1) = input_info.GEO.E(1).er;
   mur(1) = input_info.GEO.E(1).mur;
   
%    input_info.GEO_app.E(1).r = R(1);
%    input_info.GEO_app.E(1).l = L(1);
   
   lambda = c0./input_info.MODES.f_axis(end);
   
   dl = lambda./10;
   
   for i = 2:input_info.GEO.N+1
       Len_Cone = input_info.GEO.E(i).l;
       R_Cone = input_info.GEO.E(i).r;
       R_Cone_prev = input_info.GEO.E(i - 1).r;
       
       N_ele(i) = round(Len_Cone/dl);
       
       R_intermediate = linspace(R_Cone_prev, R_Cone, N_ele(i) + 1);
       
       R_intermediate2 = R_intermediate(2:end);
       
       R = [R R_intermediate2];
       
       L = [L ones(1, length(R_intermediate2)) .* dl];
       
       er = [er ones(1, length(R_intermediate2)) * input_info.GEO.E(i).er];
       mur = [mur ones(1, length(R_intermediate2)) * input_info.GEO.E(i).mur];
       
   end
   
   for k = 1:length(R)
       input_info.GEO_app.E(k).r = R(k);
       input_info.GEO_app.E(k).l = L(k);
   end
   
   plot_app = input('Want to plot the discretized version (approximate) model of the Geometry? 1 for Yes, 0 for No: [0] ');


   if isempty(plot_app)
      plot_app = 0;
   end

   if plot_app == 1
       input_info.GEO_app.type = 2;
       input_info.GEO_app.N = length(R) - 1;
       Plot_GEO(input_info.GEO_app);
   end
   
end



n = length(R);

J = n - 1; % Number of Junctions

disp('The number of cylinderical elements to solve for General Scattering Matrix: '); disp(num2str(n));

disp('Number of waveguide Junction problems to be solved: '); disp(num2str(J));

N = round(linspace(input_info.MODES.Nbn, input_info.MODES.Ntn, n));

disp('Number of modes for all waveguides in succession: '); disp(num2str(N));

F = input_info.MODES.f_axis;

%% Frequency independent inner cross product 

parfor j = 1:J
    disp('Junction');
    disp(j);
    disp('of')
    disp(J);
   
    x_til = zeros(N(j), N(j + 1));
    x_til(:, :) = Inner_prod(1:1:N(j), 1:1:N(j + 1), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j));
    X_til(j).x_til = x_til;
end

%% Frequency loop to find the GSM of the entire structure
RLRR_TE11 = zeros(1, length(F));

parfor k = 1:length(F)
    
    disp('Frequency Iteration: ');
    disp(k);

if n == 2
    
    [STT_, STR_, SRT_, SRR_] = GSM_2cyl(1:1:N(1), 1:1:N(2), F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), X_til(1).x_til);

else

[S33, S34, S43, S44] = GSM_2cyl(1:1:N(1), 1:1:N(2), F(k), R(2), R(1), er(2), mur(2), er(1), mur(1), X_til(1).x_til);
[S11, S12, S21, S22] = GSM_2cyl(1:1:N(2), 1:1:N(3), F(k), R(3), R(2), er(3), mur(3), er(2), mur(2), X_til(2).x_til);
Sl = SL(R(2), F(k), 1:1:N(2), L(2));

  
[STT_, STR_, SRT_, SRR_] = cascade_3(1:1:N(2), S11, S12, S21, S22, S33, S34, S43, S44, Sl);

% Use the for loop in case of more than 3 junctions (J > 3)

for j = 3:J

    % recursion 
    
    [S11, S12, S21, S22] = GSM_2cyl(1:1:N(j), 1:1:N(j + 1), F(k), R(j + 1), R(j), er(j + 1), mur(j + 1), er(j), mur(j), X_til(j).x_til);
    S33 = STT_; S34 = STR_; S43 = SRT_; S44 = SRR_;
    Sl = SL(R(j), F(k), 1:1:N(j), L(j));
    
    [STT_, STR_, SRT_, SRR_] = cascade_3(1:1:N(j), S11, S12, S21, S22, S33, S34, S43, S44, Sl);
    
end

end

slr = SL(R(1), F(k), 1:1:N(1), L(1));
slt = SL(R(end), F(k), 1:1:N(end), L(end));

STT(k, :, :) = slt * STT_ * slt; 
STR(k, :, :) = slt * STR_ * slr; 
SRT(k, :, :) = slr * SRT_ * slt; 
SRR(k, :, :) = slr * SRR_ * slr;
SRR_1 = squeeze(SRR(k, :, :));

f_base =  fc(R(1), er(1), mur(1));
Num_modes_prop  =  find(f_base < F(k));

% RLRR_TE11(k) = db(sum(sum(abs(SRR(k, 1:Num_modes_prop(end), 1:Num_modes_prop(end))).^2)))./2; % Return loss at waveguide R
RLRR_TE11(k) = db(abs(sum(SRR_1(1, 1:Num_modes_prop(end)))).^2)./2; % Return loss at waveguide R
% RLRR_TE11(k) = db(sum(abs(SRR_1(1, :))));
% RLRR_TE11(k) = db(abs(sum(SRR(k, 1, 1:Num_modes_prop(end))))).^2./2; % Return loss at waveguide R
end



end
