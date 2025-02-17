%% Mutual admittance function for modes like TE1m and TE1n for a clyndrical waveguide

% Inputs:
%       m:          - Second mode index of first TE mode TE1m
%       n:          - Second mode index of second TE mode TE1n
%       K0:         - Wavenumber of operation
%       R:          - Radius of the waveguide
%       er:         - Relative permittivity
%       mur:        - Relative permeability
%       L:          - Integration limit in multiples of K0
%Outputs:
%       Yin:        - Mutual reactance between mode TE1m and TE1n

% References:
% [1]﻿Dash, T. (2020). Computationally Efficient Conical Horn Antenna Design 
% [Delft University of Technology]. http://resolver.tudelft.nl/uuid:190e87c7-9309-470f-a821-43b7c3b8867b

function [Yin] = Yin_Circular_TE(m, n, k0, R, er, mur, L) 

c0 = 3e8;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er .* er0;
mu = mur .* mu0;

omega = c0 .* k0;

F = omega./(2 * pi);


% beta = 2 * pi * F ./ c0;

omega = 2 * pi * F;

Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;

xmn = str(m).xmn;
% m = str(m).m;
mode = str(m).mode;

xmn1 = str(n).xmn;
% m1 = str(n).m;
mode1 = str(n).mode;

beta_rhop = xmn./R;
% beta_z = -1j .* sqrt(-(beta.^2 - beta_rhop.^2));

beta_rhop1 = xmn1./R;
% beta_z1 = -1j .* sqrt(-(beta.^2 - beta_rhop1.^2));

% ZTE = 2 * pi * F * mu ./ beta_z;
% YTE = 1./ZTE;
% ZTE1 = 2 * pi * F * mu ./ beta_z1;


% if m == 0
%   deltam = 1;
% else
deltam = 0;
% end

Nup = (pi*(1+deltam)/2 .* (xmn.^2 - 1.^2) .* (besselj(1, xmn)).^2).^(-1);
Nup1 = (pi*(1+deltam)/2 .* (xmn1.^2 - 1.^2) .* (besselj(1, xmn1)).^2).^(-1);


k0 = (2 * pi * F)/c0;
%%

% Nu = linspace(0, 50.*beta_rhop1, 51);

% Nu = beta_rhop1;
% dnu = Nu(2) - Nu(1);

kz = @(Nu) -1j .* sqrt(-(k0.^2 - Nu.^2));

% for i = 1:length(Nu)

% ep = 1e-7;
ep = eps;

if m == n
    
    
I0 = @(Nu) Lommel2(0, R, beta_rhop, Nu, 0, 0);
I2 =  @(Nu) Lommel2(0, R, beta_rhop, Nu, 2, 2);
        
I0_2 = @(Nu) Lommel2(0, R, beta_rhop1, Nu, 0, 0);
I2_2 = @(Nu) Lommel2(0, R, beta_rhop1, Nu, 2, 2);

  I0_ = @(Nu) 1./2 .* R^2 .* (besselj(0, Nu .* R).^2 - ...
besselj(-1, Nu .* R) .*  besselj(1, Nu .* R)); 
     I2_ = @(Nu) 1./2 .* R^2 .* (besselj(2, Nu .* R).^2 - ...
besselj(1, Nu .* R) .*  besselj(3, Nu .* R)); 
   
I =  @(Nu) sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* pi ./ (4 .* omega .* mu .* kz(Nu)) .* ((2 .* k0.^2 - Nu.^2).*...
(I0(Nu) .* I0_2(Nu) + I2(Nu) .* I2_2(Nu)) - Nu.^2 .* (I0(Nu) .* I2_2(Nu) + I0_2(Nu) .* I2(Nu))) .* Nu ;
    
   
% I =  @(Nu) sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* pi ./ (4 .* omega .* mu .* kz(Nu)) .* ((k0.^2 - Nu.^2).*...
% (I0(Nu) .* I0_2(Nu) + I2(Nu) .* I2_2(Nu)) - Nu.^2 .* 1./2 .* (I0(Nu) .* I2_2(Nu) + I0_2(Nu) .* I2(Nu))) .* Nu ;


I_near_betarho = @(Nu) sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* pi ./ (4 .* omega .* mu .* kz(Nu)) .* ((2 .* k0.^2 - Nu.^2).*...
(I0_(Nu) .* I0_(Nu) + I2_(Nu) .* I2_(Nu)) - Nu.^2 .* (I0_(Nu) .* I2_(Nu) + I0_(Nu) .* I2_(Nu))) .* Nu ;

% I_near_betarho = @(Nu) sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* pi ./ (4 .* omega .* mu .* kz(Nu)) .* ((k0.^2 - Nu.^2).*...
% (I0_(Nu) .* I0_(Nu) + I2_(Nu) .* I2_(Nu)) - Nu.^2 .* 1./2 .* (I0_(Nu) .* I2_(Nu) + I0_(Nu) .* I2_(Nu))) .* Nu ;

    Y1 = integral(I, eps, beta_rhop-ep+1j.*ep, 'Waypoints', [ep+1j.*ep]);
    Y2 = integral(I_near_betarho, beta_rhop-ep+1j.*ep, beta_rhop+ep+1j.*ep);
    Y3 = integral(I, beta_rhop+ep+1j.*ep, L.*k0);
    Yin = Y1 + Y2 + Y3;
%     Yin = integral(I, eps, L.*k0); % 'Waypoints', wpoints2);

else
    beta_rhopmin = min(beta_rhop, beta_rhop1);
    beta_rhopmax = max(beta_rhop, beta_rhop1);
    
    
I0 = @(Nu) Lommel2(0, R, beta_rhopmin, Nu, 0, 0);
I2 =  @(Nu) Lommel2(0, R, beta_rhopmin, Nu, 2, 2);
        
I0_2 = @(Nu) Lommel2(0, R, beta_rhopmax, Nu, 0, 0);
I2_2 = @(Nu) Lommel2(0, R, beta_rhopmax, Nu, 2, 2);
   
I0_ = @(Nu) 1./2 .* R^2 .* (besselj(0, Nu .* R).^2 - ...
besselj(-1, Nu .* R) .*  besselj(1, Nu .* R)); 
I2_ = @(Nu) 1./2 .* R^2 .* (besselj(2, Nu .* R).^2 - ...
besselj(1, Nu .* R) .*  besselj(3, Nu .* R)); 

I =  @(Nu) sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* pi ./ (4 .* omega .* mu .* kz(Nu)) .* ((2 .* k0.^2 - Nu.^2).*...
(I0(Nu) .* I0_2(Nu) + I2(Nu) .* I2_2(Nu)) - Nu.^2 .* (I0(Nu) .* I2_2(Nu) + I0_2(Nu) .* I2(Nu))) .* Nu ;
    
I_near_betarhomin = @(Nu) sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* pi ./ (4 .* omega .* mu .* kz(Nu)) .* ((2 .* k0.^2 - Nu.^2).*...
(I0_(Nu) .* I0_2(Nu) + I2_(Nu) .* I2_2(Nu)) - Nu.^2 .* (I0_(Nu) .* I2_2(Nu) + I0_2(Nu) .* I2_(Nu))) .* Nu ;

I_near_betarhomax =  @(Nu) sqrt(Nup) .* sqrt(Nup1) .* beta_rhop .* beta_rhop1 .* pi ./ (4 .* omega .* mu .* kz(Nu)) .* ((2 .* k0.^2 - Nu.^2).*...
(I0(Nu) .* I0_(Nu) + I2(Nu) .* I2_(Nu)) - Nu.^2 .* (I0(Nu) .* I2_(Nu) + I0_(Nu) .* I2(Nu))) .* Nu ;

    Y1 = integral(I, eps, beta_rhopmin-ep+1j.*ep, 'Waypoints', [ep+1j.*ep]);
    Y2 = integral(I_near_betarhomin, beta_rhopmin-ep+1j.*ep, beta_rhopmin+ep+1j.*ep);
    Y3 = integral(I, beta_rhopmin+ep+1j.*ep, beta_rhopmax-ep+1j.*ep);
    Y4 = integral(I_near_betarhomax, beta_rhopmax-ep+1j.*ep, beta_rhopmax+ep+1j.*ep);
    Y5 = integral(I, beta_rhopmax+ep+1j.*ep, L.*k0);
    
    Yin = Y1 + Y2 + Y3 + Y4 + Y5;
end
%E_2 = integral(I, eps, k0+eps);

% Yin = E_1./YTE;
% 
% Gamma_K = (1 - Yin)./(1 + Yin);

end

