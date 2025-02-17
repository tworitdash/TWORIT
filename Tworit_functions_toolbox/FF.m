%% Far Field computation of a cascaded cylindrical waveguide/ conical waveguide problems.

% Inputs: 
%       input_info:           It is a structure input having the following fileds
%       plot_info.GEO -       Geometrical data of the waveguide structure
%       plot_info.MODE  -     Generally contains the number of modes the
%                             base and the top waveguide have
%       plot_info.FF  -       It has all Far field parameters including
%                             frequency axis (f_axis), elevation and azimuth
%                             intervals (th, ph), the general scattering matrix for
%                             this kind of geometry (GSM), the normalization
%                             procedure for far fields (norm), the excitation
%                             vector at the base waveguide (ar), the excitation
%                             modes (ex_modes), the excitation mode indices
%                             (mff)
% Outputs: 
%        Eth:                 Elevation component of the electric far field
%        Eph:                 Azimuth component of the electric far field
%        E:                   Total electric far field
%        Eco:                 Co-polar electric far field
%        Exp:                 Cross-Polar electric far field
%                             Dimensions of all the field outputs:
%                             DIM1: Frequencies
%                             DIM2: Azimuthal angle
%                             DIM3: Elevation angle


function [Eth, Eph, Eco, Exp, E, Gamma, Dm, R11, Rap] = FF(input_info)
    c0 = 3e8; % Speed of light
    modes = input_info.FF.ex_mode;
    ar = zeros(size(input_info.FF.GSM.SRR, 2), 1);
    ar(input_info.FF.mff) = input_info.FF.ar;
    
    str = input_info.FF.GSM.STR;
    stt = input_info.FF.GSM.STT;
    
    modes = load('Xmn.mat');
    modes = modes.Xmn;
    
    Ntn = input_info.MODES.Ntn;
    modest = modes(1:Ntn);
    
    norm = input_info.FF.norm;
    ap = input_info.FF.ap;
    
    rho_ = linspace(eps, input_info.GEO.E(end).r, 100);
    phi_ = linspace(eps, 2*pi, 360);
    
    drho = rho_(2) - rho_(1); dphi = phi_(2) - phi_(1);
    
    [rho, phi] = meshgrid(rho_, phi_);
    
    Erhof = zeros(length(input_info.FF.f_axis), size(rho, 1), size(rho, 2));
    Ephif = zeros(length(input_info.FF.f_axis), size(rho, 1), size(rho, 2));
    
    %% Aperture reflection [if it is set from the options]
    
    at = zeros(length(input_info.FF.f_axis), size(input_info.FF.GSM.STT, 2));
    
    if ap == 1
        % if the option to consider the aperture reflection is set (ap ==
        % 1), this sub-routine is used 
        
        [Gamma, Dm, R11, Rap, I] = AperRef(input_info.GEO.E(end).r, input_info.FF.f_axis, Ntn, modest, input_info.GEO.E(end).er, input_info.GEO.E(end).mur);
    
        at(I, :) = Gamma;
    end
    
    %%
    
    for f = 1:length(input_info.FF.f_axis)
        bt = squeeze(stt(f, :, :)) * at(:, f) + squeeze(str(f, :, :)) * ar;
        [Erho, Ephi] = EF(Ntn, rho, phi, modest, norm, input_info.FF.f_axis(f), input_info.GEO.E(end).r, 0, input_info.GEO.E(end).er, input_info.GEO.E(end).mur);
        for i = 1:Ntn
            Erhof(f, :, :) = squeeze(Erhof(f, :, :)) + bt(i) .* squeeze(Erho(i, :, :));
            Ephif(f, :, :) = squeeze(Ephif(f, :, :)) + bt(i) .* squeeze(Ephi(i, :, :));
        end
    end
    
    for f = 1:length(input_info.FF.f_axis)
        Ex(f, :, :) = cos(phi) .* squeeze(Erhof(f, :, :)) - sin(phi) .* squeeze(Ephif(f, :, :));
        Ey(f, :, :) = sin(phi) .* squeeze(Erhof(f, :, :)) + cos(phi) .* squeeze(Ephif(f, :, :));
    end
    
    
    th_ = input_info.FF.th .* pi/180;
    ph_ = input_info.FF.ph .* pi/180;
    
    [th, ph] = meshgrid(th_, ph_);
    
    for f = 1:length(input_info.FF.f_axis)
        k0 = 2.*pi.* input_info.FF.f_axis(f)./c0;
        for t = 1:length(th_)
            for p = 1:length(ph_)
                Ix_tp = squeeze(Ex(f, :, :)) .* exp(-1j .* k0 .* rho .* sin(th_(t)) .* cos(ph_(p) - phi)) .* rho .* drho .* dphi;
                Ix(p, t) = sum(sum(Ix_tp));
                
                Iy_tp = squeeze(Ey(f, :, :)) .* exp(-1j .* k0 .* rho .* sin(th_(t)) .* cos(ph_(p) - phi)) .* rho .* drho .* dphi;
                Iy(p, t) = sum(sum(Iy_tp));
                
            end
            
        end
        Eth(f, :, :) = Ix .* cos(ph) + Iy .* sin(ph);
        Eph(f, :, :) = cos(th) .* (-Ix .* sin(ph) + Iy .* cos(ph));
        Eco(f, :, :) = sin(ph) .* squeeze(Eth(f, :, :)) + cos(ph) .* squeeze(Eph(f, :, :));
        Exp(f, :, :) = cos(ph) .* squeeze(Eth(f, :, :)) - sin(ph) .* squeeze(Eph(f, :, :));
        
        E(f, :, :) = sqrt(abs(squeeze(Eth(f, :, :))).^2 + abs(squeeze(Eph(f, :, :))).^2);
    end
end