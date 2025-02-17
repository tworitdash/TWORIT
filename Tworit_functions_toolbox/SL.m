function [SL] = SL(rr, F, Nr, L)
%% SL Returns the phase terms of a bridging waveguide with radius rr and length L
% Inputs:
%        rr                 - Radius of the cylindrical waveguide
%        F                  - Frequency of Operation
%        Nr                 - Vector of mode numbers in the second
%                             waveguide (bridging waveguide between 1st and 3rd waveguide)
%                             Example: linspace(1, 10, 10)
% Outputs:      
%        SL                 - Phase term with a diagonal matrix of size 
%                             length(Nr)-by-length(Nr)
%                             containing the phase of the EM waves
%                             encountered in the second waveguide (bridging
%                             waveguide) sl = diag(e^{-1j \beta_{zn} l})
%                             where $\beta_{zn}$ is the propagation
%                             wavenumber of the nth mode propagating in
%                             waveguide no 2 and l is the length of the
%                             waveguide no 2.                                                  
% References: 
% [1]﻿Dash, T. (2020). Computationally Efficient Conical Horn Antenna Design 
% [Delft University of Technology]. http://resolver.tudelft.nl/uuid:190e87c7-9309-470f-a821-43b7c3b8867b
% [2] Tak Sum Chu and Itoh, T. (1986). Generalized Scattering Matrix Method for Analysis of Cascaded and Offset Microstrip Step Discontinuities. 
% IEEE Transactions on Microwave Theory and Techniques, 34(2):280–284.
%%
    c0 = 3e8;
    
    
    Str = load('Xmn.mat');
    Xmn = Str.Xmn;
   
    SL = zeros(length(Nr), length(Nr));
    
    for i = 1:length(Nr)
        beta = (2 * pi * F) ./ c0;
        beta_rho = (Xmn(i).xmn)./rr;
        
        beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
        
        SL(i, i) = exp(-1j .* beta_z .* L);
    end
end