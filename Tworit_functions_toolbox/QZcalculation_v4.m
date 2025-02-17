
function [Q, Z, Y] = QZcalculation_v4(N, F, r, er, mur) % , rho_, phi_, z, drho, dphi)

%% QZcalculation_v4 returns the Normalization, Impedance and Admittance matrix of a cylindrical waveguide with modes in rows and columns:
% Inputs:
%       N -   Vector of mode numbers of the cylindrical waveguide
%       F -   Frequency of operation
%       r -   Radius of the cylindrical waveguide
%       er -  Relative permittivity 
%       mur - Relative permeability
% Outputs:
%       Q -   Normalization constants in a diagonal matrix: length(N)-by-length(N)
%       Z -   Impedance diagnoal matrix                   : length(N)-by-length(N)
%       Y -   Admittance diagonal matrix                  : length(N)-by-length(N)
% References:
% [1]﻿Dash, T. (2020). Computationally Efficient Conical Horn Antenna Design 
% [Delft University of Technology]. http://resolver.tudelft.nl/uuid:190e87c7-9309-470f-a821-43b7c3b8867b
%%

c0 = 3e8;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium
mu = mu0 * mur;


% Q = zeros(size(F, 2), N(end), N(end));
% Q = zeros(length(N), length(N));
% Z = zeros(length(N), length(N));
% Y = zeros(length(N), length(N));


% xmn_ = zeros(length(N));
% beta_z_ = zeros(N);

% K = zeros(length(N), length(N));


%Y = zeros(m(end), m(end));

omega = 2 * pi * F;

lamb = c0./F; % wavelength
 
beta = 2 * pi ./ lamb;


% f =  fc(r, er, mur);
% N_i  =  find(f < F);

Str = load('Xmn.mat');
Xmn = Str.Xmn;

x = vertcat(Xmn(1:1:length(N)).xmn);
beta_rho = x./r;
beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));

Q = beta_z./abs(beta_z);

for i = 1:length(N)

   
      
    mode = Xmn(N(i)).mode;
    
   if mode == "TE"
        Z_i = 2 * pi * F * mu./ beta_z(i);
        Y_i = 1./Z_i;
%         Q_i = 1./Z_i;
%         K_i = beta_z ./ (omega .* mu .* epsilon^2);
    elseif mode == "TM"
        Z_i = beta_z(i) ./ (2 * pi * F .* epsilon);
        Y_i = 1./Z_i;
%         Q_i = 1./Z_i;
%         K_i = beta_z ./ (omega .* mu.^2 .* epsilon);
   end

    
    Z(i) = Z_i;
    Y(i) = Y_i;
%     Q(i) = Q_i;
%     K(i) = K_i;
    
end

Q = diag(Q);
Z = diag(Z);
Y = diag(Y);
% K = diag(K);


end