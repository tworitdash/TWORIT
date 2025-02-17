function [STT, STR, SRT, SRR] = cascade_3(N, S11, S12, S21, S22, S33, S34, S43, S44, Sl)

%% cascade_3 Returns the S parameters of a 3 cylindrical waveguide structre 
% Inputs:
%        N -                  Number of modes in the second waveguide (bridging waveguide between 1st and 3rd waveguide)
%        S11, S12, S21, S22 - S parameters from the General Scattering
%                             Matrix (GSM) of Waveguide no 3 and waveguide
%                             no 2 - Check GSM_2cyl function for more
%                             information about how to find these matrices
%        S33, S34, S43, S44 - S parameters from the General Scattering
%                             Matrix (GSM) of Waveguide no 2 and waveguide no 1
%                             Check GSM_2cyl function for more
%                             information about how to find these matrices
%        Sl                 - Phase term with N by N diagonal matrix
%                             containing the phase of the EM waves
%                             encountered in the second waveguide (bridging
%                             waveguide) sl = diag(e^{-1j \beta_{zn} l})
%                             where $\beta_{zn}$ is the propagation
%                             wavenumber of the nth mode propagating in
%                             waveguide no 2 and l is the length of the
%                             waveguide no 2.
%                             Check SL function to find out the code for
%                             this
% Outputs:
%       [STT, STR, SRT, SRR] - Is nothing but [S11, S14, S41, S44] T is
%                              usually used for the bigger waveguide and R is used for the smaller
%                              waveguide. 
% References: 
% [1]﻿Dash, T. (2020). Computationally Efficient Conical Horn Antenna Design 
% [Delft University of Technology]. http://resolver.tudelft.nl/uuid:190e87c7-9309-470f-a821-43b7c3b8867b
% [2] Tak Sum Chu and Itoh, T. (1986). Generalized Scattering Matrix Method for Analysis of Cascaded and Offset Microstrip Step Discontinuities. 
% IEEE Transactions on Microwave Theory and Techniques, 34(2):280–284.
%%
    I = eye(length(N), length(N));
    
    U1 = inv(I - S22 * Sl * S33 * Sl);
    U2 = inv(I - S33 * Sl * S22 * Sl);

    STT = (S11 + S12 * Sl * U2 * S33 * Sl * S21);
    STR = (S12 * Sl * U2 * S34);
    SRT = (S43 * Sl * U1 * S21);
    SRR = (S44 + S43 * Sl * U1 * S22 * Sl * S34);

end