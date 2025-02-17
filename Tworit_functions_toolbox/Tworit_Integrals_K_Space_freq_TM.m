function [Dm_, yap, Gamma, y11, Gamma_11, YTM11, y_for_debug,Ymut_k,Y_rhs_k] = Tworit_Integrals_K_Space_freq_TM(r, N, k0, L, er, mur, n_orig)

c0 = 3e8;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

omega = c0.*k0;

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

Gamma_11(k) = (1 - 1./y11(k))./(1 + 1./y11(k));

Gamma(k) = (1 - 1./yap(k))./(1 + 1./yap(k));

% Gamma_11(k) = (1 + y11(k))./(y11(k) - 1);
% 
% Gamma(k) = (1 + yap(k))./(yap(k) - 1);

% Gamma_11(k) = (1 - y11(k))./(1 + y11(k));
% % 
% Gamma(k) = (1 - yap(k))./(1 + yap(k));
% 


end
end