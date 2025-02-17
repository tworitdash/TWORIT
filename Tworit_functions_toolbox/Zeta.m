function [Xx] = Zeta(x, beta_rhoNu, beta_rhoMu, Nu, Mu)

 Xx = (x./(beta_rhoNu.^2 - beta_rhoMu.^2)) .* (-beta_rhoNu .* besselj(Nu - 1, beta_rhoNu .* x)...
 .* besselj(Mu, beta_rhoMu .* x) + beta_rhoMu .*  besselj(Nu, beta_rhoNu .* x)...
 .* besselj(Mu - 1, beta_rhoMu .* x));
    
end