function dJmdZ = besselj_der(m, z)

%% Calculate the Derivative of the Bessel function 
% besselj_der(m, z) returns The derivative of the Bessel function of order
% 'm' and argument z

dJmdZ= (besselj(m - 1, z) - besselj(m + 1, z))./2;

end