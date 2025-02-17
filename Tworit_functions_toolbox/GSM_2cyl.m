function [Spp, Spr, Srp, Srr] = GSM_2cyl(Nr, Np, F, rp, rr, erp, murp, err, murr, X_)

%% GSM_2cyl function returns the General Scattering Matrix (GSM) of the 2 cylinder
% junction problem 
% Inputs:
%       Nr -   Mode numbers of the smaller waveguide example: Nr =
%              linspace(1, 10, 10); The modes corresponding to serial
%              numbers can be found in MODES.csv
%       Np -   Mode numbers of the bigger waveguide example:  Np =
%              linspace(1, 10, 40); The modes corresponding to serial
%              numbers can be found in MODES.csv
%       F -    Frequency of operation
%       rp -   Radius of the bigger waveguide in [m]
%       rr -   Radius of the smaller waveguide in [m]
%       erp -  Relative permittivity of the bigger waveguide
%       murp - Relative permeability of the bigger waveguide
%       err -  Relative permittivity of the smaller waveguide
%       murr - Relative permeability of the smaller waveguide
%       X_ -   Is a matrix of length length(Nr) by length(Np) - Known as the
%              Inner cross product between the modes of the two waveguides
%              - Check the function Inner_prod for more details about the
%              formation of this matrix
% Outputs:
%       Spp - Reflection coefficients of the modes originatig from the
%             bigger waveguide - size - length(Np) by length(Np)
%       Spr - Transmission coefficients of the modes originatig from the
%             bigger waveguide towards the modes of the smaller waveguide
%             - size - length(Np) by length(Nr)
%       Srp - Transmission coefficients of the modes originatig from the
%             smaller waveguide towards the modes of the bigger waveguide
%             - size - length(Nr) by length(Np)
%       Spp - Reflection coefficients of the modes originatig from the
%             smaller waveguide - size - length(Nr) by length(Nr)
% Reference: 
% [1]﻿Dash, T. (2020). Computationally Efficient Conical Horn Antenna Design 
% [Delft University of Technology]. http://resolver.tudelft.nl/uuid:190e87c7-9309-470f-a821-43b7c3b8867b


Spp = zeros(length(Np), length(Np));
Spr = zeros(length(Np), length(Nr));
Srp = zeros(length(Nr), length(Np));
Srr = zeros(length(Nr), length(Nr));

%% Frequency independent Modular inner cross product between the two wavegudies

X_til = zeros(length(Nr), length(Np));

 for p = 1:length(Np)
     for r = 1:length(Nr)
           X_til(r, p) = X_(r, p);
     end
 end
    
    

%% Wavwguide p - Calculation of Normalization, Impedance, Admittance

drho = rp/100;
dphi = pi/180;

[rho_, phi_] = meshgrid(eps:drho:rp, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zp = 0; 

[Qp, Zp, Yp] = QZcalculation_v4(Np, F, rp, erp, murp); %, rho_, phi_, zp, drho, dphi);


%% Wavwguide r - Calculation of Normalization, Impedance, Admittance

drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

[Qr, Zr, Yr] = QZcalculation_v4(Nr, F, rr, err, murr); % , rhor_, phir_, zr, drho, dphi);

%% Calculation of GSM

Ip = eye(length(Np), length(Np));
Ir = eye(length(Nr), length(Nr));

X = sqrt(Qr * Zr) * X_til * sqrt(Yp * Qp); % modular inner cross product. Takes the dimension of Np \times Nr

F_ = 2 * inv(Qr + X * inv(Qp) * X.');

Spp = inv(Qp) * X.' * F_ * X - Ip;
Spr = inv(Qp) * X.' * F_ * Qr;
Srp = Spr.';
Srr = F_ * Qr - Ir;


end


