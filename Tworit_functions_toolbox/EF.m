%% Electric near field computation of cylindrical waveguide

% Inputs: 
%        Ntn:                 Number of waveguide modes needed
%        rho:                 radial domain
%        phi:                 azimuthal domain
%        norm:                normalization procedure 1 is power
%                             normalization, 2 is FEKO style normalization
%        f:                   Frequency of operation
%        modest:              Details of all modes of which Ntn number of
%                             modes are used for computation
%                             It is a structural variable having variables
%                             such as, type of mode (TE/ TM), the mode
%                             numbers m and n and the polarization, 0 or 90
%        r:                   radius of the waveguide: The final value of
%                             rho
%        z:                   Z axis length where the field has to be
%                             computed
%        er:                  Relative permittivity of the last element in
%                             the Geometry
%        mur:                 Relative permeability of the last element in
%                             the Geometry
% Outputs: 
%        Erho:                Radial component of the electric near field
%                             It has the following dimensions:
%                             DIM1: Mode
%                             DIM2: radius
%                             DIM3: azimuthal angle
%        Ephi:                Azimuthal component of the electric near field
%                             It has the following dimensions:
%                             DIM1: Mode
%                             DIM2: Radius
%                             DIM3: Azimuthal angle



function [Erho, Ephi] = EF(Ntn, rho, phi, modest, norm, f, r, z, er, mur)
    for i = 1:Ntn
        [Erho(i, :, :), Ephi(i, :, :)] = Erhophi(modest(i), rho, phi, norm, f, r, z, er, mur);
    end
end