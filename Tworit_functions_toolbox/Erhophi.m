%% Electric near field components of cylindrical waveguide

% Inputs: 
%       
%        rho:                 radial domain
%        phi:                 azimuthal domain
%        norm:                normalization procedure 1 is power
%                             normalization, 2 is FEKO style normalization
%        f:                   Frequency of operation
%        modest:              Details of the mode
%                             It is a structural variable having variables
%                             such as, type of mode (TE/ TM), the mode
%                             numbers m and n and the polarization, 0 or 90
%        r:                   radius of the waveguide: The final value of
%                             rho
%        z:                   Z axis length where the field has to be
%                             computed
%        er:                  Relative permittivity of the waveguide
%        mur:                 Relative permeability of the waveguide             
% Outputs: 
%        Erho:                Radial component of the electric near field
%                             It has the following dimensions:
%                          
%                             DIM1: radius
%                             DIM2: azimuthal angle
%        Ephi:                Azimuthal component of the electric near field
%                             It has the following dimensions:
%                             
%                             DIM1: Radius
%                             DIM2: Azimuthal angle




function [Erho, Ephi] = Erhophi(modest, rho, phi, norm, f, r, z, er, mur)
c0 = 3e8;
epsilon = er .* 8.85418782e-12; % Free space permittivity
mu = mur * 1.25663706e-6;  % Free Space Permeability
 
    if modest.mode == "TE"
            A = 1;
            if modest.pol == 0
                C = 1; 
                D = 0; 
            else
                C = 0; 
                D = 1; 
            end
            
            m = modest.m;

            omega = 2 * pi * f;
            beta = omega/c0;
            
            beta_rho = modest.xmn/r;
    
            beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
   
            ZTE = (omega .* mu ./ beta_z);
    
            Erho = -1j .* A .* m .* beta_z ./(beta_rho.^2 .* rho) .* ZTE .* besselj(m, beta_rho .* rho) .* (-C .* sin(m .* phi)...
                + D .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
            Ephi = 1j .* A .* beta_z ./ beta_rho  .* ZTE .* besselj_der(m, beta_rho .* rho) .* (C .* cos(m .* phi)...
                + D .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
            Ez = zeros(size(rho)) + eps;

            Hrho = -Ephi/ZTE;
            Hphi = Erho/ZTE;

            Hz = A .* besselj(m, beta_rho .* rho) .* cos(m .* phi);
       
    else
   
            B = 1;
            if modest.pol == 0
                C = 1; 
                D = 0; 
            else
                C = 0; 
                D = 1; 
            end
            m = modest.m;
            
            omega = 2 * pi * f;
            beta = omega/c0;
            
            beta_rho = modest.xmn/r;
            
            beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
            ZTM = beta_z./(omega .* epsilon);

            Erho = -1j .*B * beta_z./beta_rho .* besselj_der(m, beta_rho .* rho).* (C .* cos(m .* phi)...
                + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
            Ephi = -1j .* B .* (m .* beta_z ./ (beta_rho.^2 .* rho)) .* besselj(m, beta_rho .* rho) .* (- C .* sin(m .* phi)...
                + D .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
            Ez = -1j .* B .* besselj(m, beta_rho .* rho) .*  (C .* cos(m .* phi)...
                + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z) ;

            Hrho = -Ephi./ZTM;
            Hphi = Erho./ZTM;
            Hz = zeros(size(rho)) + eps;
        
    end
end