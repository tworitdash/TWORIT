function [X_til] = Inner_prod(Nr, Np, rp, rr, erp, murp, err, murr) 
%% Inner_prod function retruns the frequency independent Inner cross product of
% a two cylindrical waveguide junction problem.
% Inputs:
%       Nr -   Mode numbers of the smaller waveguide example: Nr =
%              linspace(1, 10, 10); The modes corresponding to serial
%              numbers can be found in MODES.csv
%       Np -   Mode numbers of the bigger waveguide example:  Np =
%              linspace(1, 10, 40); The modes corresponding to serial
%              numbers can be found in MODES.csv
%       rp -   Radius of the bigger waveguide in [m]
%       rr -   Radius of the smaller waveguide in [m]
%       erp -  Relative permittivity of the bigger waveguide
%       murp - Relative permeability of the bigger waveguide
%       err -  Relative permittivity of the smaller waveguide
%       murr - Relative permeability of the smaller waveguide
% Output:
%       X_til - Inner cross product between the modes of the two
%               waveguides: size - length(Nr) by length(Np) 
% Reference: 
% [1]﻿Dash, T. (2020). Computationally Efficient Conical Horn Antenna Design 
% [Delft University of Technology]. http://resolver.tudelft.nl/uuid:190e87c7-9309-470f-a821-43b7c3b8867b

%% Inner Product Calculation


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;       % Permeability in the medium
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;


Str = load('Xmn.mat');
Xmn = Str.Xmn;


drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide

X_til = zeros(length(Nr), length(Np));


for p = 1:length(Np)
      for r = 1:length(Nr)
%                 disp('Iteration:')
% %                 
%                 disp(p);
%                 disp(r);
                
                xmn_p = Xmn(p).xmn;
                xmn_r = Xmn(r).xmn;
                
                beta_rhop = xmn_p/rp;
                beta_rhor = xmn_r/rr;
                
                pm = Xmn(p).m;
                rm = Xmn(r).m;
                
                modep = Xmn(p).mode;
                moder = Xmn(r).mode;
                
                polp = Xmn(p).pol;
                polr = Xmn(r).pol;
                
                if pm == 0
                    deltapm = 1;
                else
                    deltapm = 0;
                end
                 
                if rm == 0
                    deltarm = 1;
                else
                    deltarm = 0;
                end
                
                if modep == "TE"
                    Nup = (pi*(1+deltapm)/2 .* (xmn_p.^2 - pm.^2) .* (besselj(pm, xmn_p)).^2).^(-1);
                elseif modep  == "TM"
                    Nup = (pi*(1+deltapm)/2 .* (xmn_p).^2 .* (besselj_der(pm, xmn_p)).^2).^(-1);
                end
                
                if moder  == "TE"
                    Nur = (pi*(1+deltarm)/2 .* (xmn_r.^2 - rm.^2) .* (besselj(rm, xmn_r)).^2).^(-1);
                elseif moder  == "TM"
                    Nur = (pi*(1+deltarm)/2 .* (Xmn(r).xmn).^2 .* (besselj_der(rm, xmn_r)).^2).^(-1);
                end
%                 Nup = 1;
%                 Nur = 1;
%                 
%              
                if polp == 0
                    grad_Phi_rhop = sqrt(Nup) .* cos(pm .* phir_) .* besselj_der(pm, beta_rhop .* rhor_) .* beta_rhop;
                    grad_Phi_phip = (-1./rhor_) .* sqrt(Nup) .* pm .* sin(pm .* phir_) .* besselj(pm,  beta_rhop .* rhor_);
                else
                    grad_Phi_rhop = sqrt(Nup) .* sin(pm .* phir_) .* besselj_der(pm, beta_rhop .* rhor_) .* beta_rhop;
                    grad_Phi_phip = (1./rhor_) .* sqrt(Nup) .* pm .* cos(pm .* phir_) .* besselj(pm,  beta_rhop .* rhor_);
                end
%                
                
%                
              if polr == 0  
                    grad_Phi_rhor = sqrt(Nur) .* cos(rm .* phir_) .* besselj_der(rm, beta_rhor.* rhor_) .* beta_rhor;
                    grad_Phi_phir = (-1./rhor_) .* sqrt(Nur) .* rm .* sin(rm .* phir_) .* besselj(rm,  beta_rhor .* rhor_);
              else
                    grad_Phi_rhor = sqrt(Nur) .* sin(rm .* phir_) .* besselj_der(rm, beta_rhor.* rhor_) .* beta_rhor;
                    grad_Phi_phir = (1./rhor_) .* sqrt(Nur) .* rm .* cos(rm .* phir_) .* besselj(rm,  beta_rhor .* rhor_);
              end
                
    if (modep  == "TE" && moder  == "TE") || (modep  == "TM" && moder  == "TM")
            if polp == polr
                    
                    if (pm == rm) && (pm ~= 0) && (rm ~= 0)
                            A = Lommel(0, rr, beta_rhop, beta_rhor, pm - 1, rm - 1);
                      
                            D = Lommel(0, rr, beta_rhop, beta_rhor, pm + 1, rm + 1);

                    
                            Icos = intphicos(0, 2*pi, pm, rm);
                            Isin = intphisin(0, 2*pi, pm, rm);
                            K = beta_rhop .* beta_rhor./4;
                      
                            X_til_pr =  sqrt(Nup) .* sqrt(Nur) ...
                            .* K .* (A + D) .* (Icos + Isin);

                            X_til(r, p) = X_til_pr;                        
                    else 
                        
 %                           X_til_pr = (grad_Phi_rhop .* grad_Phi_rhor +  grad_Phi_phip .* grad_Phi_phir)...
  %                      .* rhor_ .* drho .* dphi;
 %                           X_til(r, p) = sum(sum(X_til_pr));

                            X_til(r, p) = 0;
%                              X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_phip)...
%                        .* rhor_ .* drho .* dphi;
%                              X_til(r, p) = sum(sum(X_til_pr));
%                   
                    end
            else

%                            X_til_pr = (grad_Phi_rhop .* grad_Phi_rhor +  grad_Phi_phip .* grad_Phi_phir)...
%                        .* rhor_ .* drho .* dphi;
 %                           X_til(r, p) = sum(sum(X_til_pr)); 
                            X_til(r, p) = 0;
            end
                    
     elseif (modep == "TE" && moder == "TM")
%                    X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_phip)...
%                      .* rhor_ .* drho .* dphi;
%                   X_til(r, p) = sum(sum(X_til_pr));
                     X_til(r, p) = 0;
                    
     elseif (modep == "TM" && moder == "TE")
                if polr == polp
                    X_til(r, p) = 0;
                end
                if (pm == rm) && (polr ~= polp)
                    X_til(r, p) = - pm .* pi .* sqrt(Nup) .* sqrt(Nur) .* (besselj(pm, beta_rhop .* rr) .* besselj(rm, beta_rhor .* rr) - besselj(pm, 0) .* besselj(rm, 0));
                else
                    X_til(r, p) = 0;
%                    X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_phip)...
%                        .* rhor_ .* drho .* dphi;
%                    X_til(r, p) = sum(sum(X_til_pr));
                end    
     end             
                    
                    
                    
        
      end
      
end
    

% csvwrite('TM_TE_Inner_P', X_til); 
% save('Inner_P_analytical_V3_ratio_1', 'X_til');