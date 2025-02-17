function [Dm_, yap, Gamma, y11, Gamma_11, YTE11, y_for_debug] = Tworit_Integrals_K_Space_freq(r, N, k0, L, er, mur, n_orig)

c0 = 3e8;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilon = er * er0;
mu = mur .* mu0;

omega = c0.*k0;

Str = load('Xmn_azimuthal_inc_TE.mat');

str = Str.xmn_TE;

for k = 1:length(k0)
    
    R = r;

for i = 1:length(N)
    
    xmn(i) = str(N(i)).xmn;
    beta_rho(i) = xmn(i)./R;
    M(i) = str(N(i)).m;
    beta_z(i) = -1j .* sqrt(-(k0(k).^2 - (beta_rho(i)).^2));
    YTE(i) = beta_z(i)./(omega(k) .* mu);
    ZTE(i) = 1./YTE(i);
    
end

for l = 1:length(N)
    for p = 1:length(N)
        if l == p
            Ymut(l, p) = (Yin_Circular(N(l), N(p), k0(k), R, er, mur, L) + YTE(p));
        else
            Ymut(l, p) = Yin_Circular(N(l), N(p), k0(k), R, er, mur, L);
        end
    end
end

for l = 1:length(N)
    Y_rhs(l) = Yin_Circular(n_orig, N(l), k0(k), R, er, mur, L);
end

Dm_(k, :) = Ymut\(-Y_rhs.');

% Dm(k, :) = Ymut\(-Y_rhs.');


% 
% xmn11 = str(1).xmn;
xmn11 = str(n_orig).xmn;
beta_rho11 = xmn11./R;

beta_z11 = -1j .* sqrt(-(k0(k).^2 - beta_rho11.^2));
YTE11(k) = beta_z11./(omega(k) .* mu);
ZTE11 = 1./YTE11(k);

% Y11 = Yin_Circular(1, 1, k0, R, er, mur, L);
Y11 = Yin_Circular(n_orig, n_orig, k0(k), R, er, mur, L);

yap(k) = (Y11 + Dm_(k, :) * Y_rhs.')./YTE11(k);
y11(k) = Y11./YTE11(k);

y_for_debug(k) = Y11;

Gamma_11(k) = (1 - y11(k))./(1 + y11(k));

Gamma(k) = (1 - yap(k))./(1 + yap(k));

end
end