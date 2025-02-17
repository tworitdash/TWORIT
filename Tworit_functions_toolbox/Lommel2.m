%% Lommel function returns the definite Lommel Integrals.
%  Input: 
%       a - Lower limit of the integration
%       b - Upper limit of the integration 
%       beta_rhoNu - Multiplicative argument of the first Bessel function 
%       beta_rhoMu - Multiplicative argument of the second Bessel function 
%       Nu - Order of the first Bessel function
%       Mu - Order of the second Bessel function NOTE: Mu and Nu are used 
%            as the same in the context of the tool 
%  Output: 
%       X : $$ \int_{a}^{b} J_{\nu}(beta_{\nu} \rho) J_{\mu}(beta_{\mu} \rho)
%           \rho d\rho $$ 
%  NOTE: 
%       Mu ~= Nu: The result corresponding to this condition is not the correct result.
%       This condition never appears in this tool, that is why it is set to 0.
% References:
% [1]﻿Dash, T. (2020). Computationally Efficient Conical Horn Antenna Design 
% [Delft University of Technology]. http://resolver.tudelft.nl/uuid:190e87c7-9309-470f-a821-43b7c3b8867b
% [2] Abramowitz, M. (1972). Handbook of Mathematical Functions with Formulas, 
% Graphs and Mathe- matical Tables. National Bureau of Standards. Applied Mathematics Series, 55 edition.
% [3] Bowman, F. (2018). Introduction to Bessel functions. Dover.



function [X] = Lommel2(a, b, beta_rhoNu, beta_rhoMu, Nu, Mu)

        Xa = Zeta(a, beta_rhoNu, beta_rhoMu, Nu, Mu);
        Xb = Zeta(b, beta_rhoNu, beta_rhoMu, Nu, Mu);
    
        X = Xb - Xa;

end