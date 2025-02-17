%% Zeros of Bessel's function (For TM) and zeros of the derivative of the Bessel's function (For TE) Calculation
clear;

m = 0:1:100;
% m = 1:1:100;
% m = 1;


for j = 1:length(m)

disp(j);

xm = linspace(0.1, 1000, 10000);

Jm = @(z) besselj(m(j), z);  % Bessel's function

dz = 1e-5;

Jm_der = @(z) (besselj(m(j), z + dz) - besselj(m(j), z - dz))./(2 * dz); % Derivative of the Bessel's function

ym_TE = Jm_der(xm); % Values of the Bessel's function at the positions defined by xm
ym_TM = Jm(xm);

if m(j) == 0
    chsign_TE = [find(diff(sign(ym_TE)))];  % Detection of the sign changes 
    chsign_TM = [find(diff(sign(ym_TM)))];
else
    chsign_TE = [find(diff(sign(ym_TE))) find(diff(sign(ym_TE)))];  % Detection of the sign changes 
    chsign_TM = [find(diff(sign(ym_TM))) find(diff(sign(ym_TM)))];
end

J_TE = 1:size(chsign_TE, 2); % number of sign changes
J_TM = 1:size(chsign_TM, 2);

if m(j) == 0
    n_TE = [J_TE(1:length(J_TE))];
    n_TM = [J_TM(1:length(J_TM))];
else
    n_TE = [J_TE(1:length(J_TE)/2) J_TE(1:length(J_TE)/2)];
    n_TM = [J_TM(1:length(J_TM)/2) J_TM(1:length(J_TM)/2)];
end

if m(j) == 0
    pol_TE = [zeros(1, length(J_TE))];
    pol_TM = [zeros(1, length(J_TM))];
else
    pol_TE = [zeros(1, length(J_TE)/2) ones(1, length(J_TE)/2) .* pi/2];
    pol_TM = [zeros(1, length(J_TM)/2) ones(1, length(J_TM)/2) .* pi/2];
end

TE_modes_(j+1) = length(J_TE);
TM_modes_(j+1) = length(J_TM);

TE_modes(j) = sum(TE_modes_(1:j));
TM_modes(j) = sum(TM_modes_(1:j));

for i = 1:TE_modes_(j+1)
    
        xmn_TE(TE_modes(j) + i).xmn = fzero(Jm_der, xm(chsign_TE(i)));  % finding the roots near the points where derivative of Jm changes sign
        xmn_TE(TE_modes(j) + i).m = m(j);
        xmn_TE(TE_modes(j) + i).n = n_TE(i);
        xmn_TE(TE_modes(j) + i).mode = "TE";
        xmn_TE(TE_modes(j) + i).pol = pol_TE(i);
        disp(i);
 
end

for k = 1:TM_modes_(j+1)
        xmn_TM(TM_modes(j) + k).xmn = fzero(Jm, xm(chsign_TM(k))); % finding the roots near the points where Jm changes sign
        xmn_TM(TM_modes(j) + k).m = m(j);
        xmn_TM(TM_modes(j) + k).n = n_TM(k);
        xmn_TM(TM_modes(j) + k).mode = "TM";
        xmn_TM(TM_modes(j) + k).pol = pol_TM(k);
        disp(k);
   
end

end

Xmn = [xmn_TE xmn_TM];


% 
[x,idx]=sort([Xmn.xmn]);
% 
Xmn = Xmn(idx);

for i = 1:length(Xmn)
    Xmn(i).sl_no = i;
end
% 
% [x_TE, idx_TE] = sort([xmn_TE.n]);
% xmn_TE = xmn_TE(idx_TE);
% 
% [x_TE, idx_TE] = sort([xmn_TE.m]);
% xmn_TE = xmn_TE(idx_TE);
% 
% 
% [x,idx]=sort([Xmn.n]);
% Xmn = Xmn(idx);
% 
% [x,idx]=sort([Xmn.m]);
% Xmn = Xmn(idx);


% [x_TM, idx_TM] = sort([xmn_TM.n]);
% xmn_TM = xmn_TM(idx_TM);
% 
% [x_TM, idx_TM] = sort([xmn_TM.m]);
% xmn_TM = xmn_TM(idx_TM);


% [x,idx]=sort([Xmn.n]);
% Xmn = Xmn(idx);
% 
% [x,idx]=sort([Xmn.m]);
% Xmn = Xmn(idx);


save('Xmn', 'Xmn');
% save('Xmn_azimuthal_inc_TE', 'xmn_TE');
% save('Xmn_azimuthal_inc_TM', 'xmn_TM');
