%% Aperture reflection function for cylindrical waveguides with modes 
%% This function is for research purposes
%% having the first order azimuthal variation like TE11, TM11, TE12


% Inputs: 
%       r               -     The radius of the waveguide
%       f_axis          -     Frequency axis
%       N               -     Number of modes to be considered
%                             
%       modes           -     Details of the N modes. It is a structural
%                             variable. For details of the fields of this
%                             structre, check Xmn.mat file in the parent
%                             directory.
%       er              -     Relative permittivity of the waveguide: This
%                             software works best for free space [er = 1]
%       mur             -     Relative permeability of the waveguide: This
%                             software works best for free space [mur = 1]
% Outputs: 
%       Gamma:          -     Reflection coeffcient at the aperture having
%                             the following dimensions:
%                             DIM1: frequency axis
%                             DIM2: Number of modes having 1 as the first
%                             index like TE11, TM11, TE12 and so on that
%                             comes under N modes
%       Dm:             -     For research purposes: Higher order mode
%                             coefficient at the boundary. 
%                             DIM1: frequency axis
%                             DIM2: Same as Gamma, number of modes having 1
%                             as first index as TE11, TM11, TE12 and so on
%                             that comes under N
%                             DIM3: 3, as we have taken 3 higher order mode
%                             exctations locally at the boundary 
%       R11:            -     Mutual reactance of each mode with itself.
%                             Dimensions same as Gamma
%       Rap:                  Aperture reactance considering all higher
%                             order modes into account
%       I:              -     Mode indices in the variable modes which are
%                             considered for aperture reflection
% NOTE: Reactance here can mean admittance for TE modes and impedance for
% TM modes


function [Gamma, Dm, R11, Rap, I] = AperRef(r, f_axis, N, modes, er, mur)
    k = 0;
    for i = 1:N
        if modes(i).m == 1
            k = k + 1;
            I(i) = i;
            b(k).mode = modes(i).mode;
            b(k).m = modes(i).m;
            b(k).n = modes(i).n;
            b(k).xmn = modes(i).xmn;
            b(k).pol = modes(i).pol;
            b(k).N = b(k).n+1:b(k).n+3;
        end
    end
    
    for l = 1:k
        if b(l).mode == 'TE'
            [Gamma(l, :), Dm(l, :, :), R11(l, :), Rap(l, :)] = Tworit_Gamma_OpenTE(r, f_axis, b(l), er, mur, 50);
        else
            [Gamma(l, :), Dm(l, :, :), R11(l, :), Rap(l, :)] = Tworit_Gamma_OpenTM(r, f_axis, b(l), er, mur, 50);
        end
    end
    
end