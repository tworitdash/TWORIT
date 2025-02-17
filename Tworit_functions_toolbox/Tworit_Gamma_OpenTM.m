%% Aperture reflection function for TM1n modes

% Inputs: 
%       r               -     The radius of the waveguide
%       f               -     Frequency axis
%                             
%       b               -     mode details. It is a structure variable. It
%                             has the same fields as Xmn.mat in the parent
%                             directory and one extra entry with indices
%                             for higher order modes in TE withe same
%                             azimuthal variation. For example, for TM11,
%                             the extry is b.N = [2 3 4] to consider TM12,
%                             TM13 and TM14
%                             
%       er              -     Relative permittivity of the waveguide: This
%                             software works best for free space [er = 1]
%       mur             -     Relative permeability of the waveguide: This
%                             software works best for free space [mur = 1]
% Outputs: 
%       Gamma:          -     Reflection coeffcient at the aperture having
%                             the following dimensions:
%                             DIM1: frequency axis
%                             
%       Dm_:             -    For research purposes: Higher order mode
%                             coefficient at the boundary. 
%                             DIM1: frequency axis                            
%                             DIM2: 3, as we have taken 3 higher order mode
%                             exctations locally at the boundary 
%       z11:            -     Mutual reactance of first mode with itself.
%                             Dimensions same as Gamma
%       zap:                  Aperture reactance considering all higher
%                             order modes into account

% References:
% [1]﻿Dash, T. (2020). Computationally Efficient Conical Horn Antenna Design 
% [Delft University of Technology]. http://resolver.tudelft.nl/uuid:190e87c7-9309-470f-a821-43b7c3b8867b


function [Gamma, Dm_, z11, zap] = Tworit_Gamma_OpenTM(r, f, b, er, mur, L)
c0 = 3e8;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

omega = 2*pi*f;
% omega = c0.*k0;
k0 = omega/c0;

N = b.N;
n_orig = b.n;

Str = load('Xmn_azimuthal_inc_TM.mat');

str = Str.xmn_TM;

for k = 1:length(k0)
    
    R = r;

for i = 1:length(N)
    
    xmn(i) = str(N(i)).xmn;
    beta_rho(i) = xmn(i)./R;
    M(i) = str(N(i)).m;
    beta_z(i) = -1j .* sqrt(-(k0(k).^2 - (beta_rho(i)).^2));
    ZTM(i) = beta_z(i)./(omega(k) .* epsilon);
    YTM(i) = 1./ZTM(i);
    
end

for l = 1:length(N)
    for p = 1:length(N)
        if l == p
            Ymut(l, p) = (Yin_Circular_TM(N(l), N(p), k0(k), R, er, mur, L) + ZTM(p));
            Ymut_k(k, l, p) = Ymut(l, p);
        else
            Ymut(l, p) = Yin_Circular_TM(N(l), N(p), k0(k), R, er, mur, L);
            Ymut_k(k, l, p) = Ymut(l, p);
        end
    end
end

for l = 1:length(N)
    Y_rhs(l) = Yin_Circular_TM(n_orig, N(l), k0(k), R, er, mur, L);
    Y_rhs_k(k, l) = Y_rhs(l);
end

Dm_(k, :) = Ymut\(-Y_rhs.');

% Dm(k, :) = Ymut\(-Y_rhs.');


% 
% xmn11 = str(1).xmn;
xmn11 = str(n_orig).xmn;
beta_rho11 = xmn11./R;


beta_z11 = -1j .* sqrt(-(k0(k).^2 - beta_rho11.^2));
ZTM11(k) = beta_z11./(omega(k) .* epsilon);
YTM11(k) = 1./ZTM11(k);

% Y11 = Yin_Circular(1, 1, k0, R, er, mur, L);
Y11 = Yin_Circular_TM(n_orig, n_orig, k0(k), R, er, mur, L);

yap(k) = (Y11 + Dm_(k, :) * Y_rhs.')./ZTM11(k);
y11(k) = Y11./ZTM11(k);

y_for_debug(k) = Y11;

z11(k) = 1./y11(k); zap(k) = 1./yap(k);

Gamma_11(k) = (1 - 1./y11(k))./(1 + 1./y11(k));

Gamma(k) = (1 - 1./yap(k))./(1 + 1./yap(k));

% Gamma_11(k) = (1 + y11(k))./(y11(k) - 1);
% 
% Gamma(k) = (1 + yap(k))./(yap(k) - 1);

% Gamma_11(k) = (1 - y11(k))./(1 + y11(k));
% % 
% Gamma(k) = (1 - yap(k))./(1 + yap(k));
end