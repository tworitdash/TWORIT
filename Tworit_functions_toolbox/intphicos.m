function [I2rho] = intphicos(a, b, pm, rm)

%% intphicos function returns the value of the definite integral 
% $$ \int_{a}^{b} cos^(pm \phi) cos(rm \phi) d\phi$$
% Inputs: 
%       a -  Lower limit of the integral
%       b -  Upper limit of the integral
%       pm - Multiplicative argument inside first cosine function 
%       pm - Multiplicative argument inside second cosine function
% Outputs:
%       I2rho : $$ \int_{a}_{b} cos^(pm \phi) cos(rm \phi) d\phi$$

%% 
if pm == rm
    I2rho = (1./2) .* ( (b - a) + (1./(2.*pm)) .* (sin(2.*pm.*b) - sin(2.*pm.*a)));
else

I2rho = (1./2) .* ( (1./(pm - rm) .* ( sin(b.*(pm - rm)) - sin(a .*(pm - rm)))) + ...
(1./(pm + rm) .* ( sin(b.*(pm + rm)) - sin(a .*(pm + rm)))) );
end

end