function [min_value,rho0,factor] = analyticerrorestimate(prefactor,lambda,phi,K,rho0,expfactor)
% Approximately minimizes (over rho>1) the quantity
% prefactor(rho) * sup_{z in E_rho} |exp(lambda*(z+1))| * sup_{z in E_rho} |phi(rho)|,
% where E_rho is the Berstein ellipse of size rho 
% and phi is represented by its Chebyshev coefficients.
% K is the leading order decay of the prefactor (prefactor(rho) ~ rho^(-K),
% which is used to get an heuristic upper bound for rho.
% In case a value of rho0 is already given in the last (optional) parameter,
% the optimization over rho is skipped, and the input value used
% Also returns the phi-independent (but rho-dependent) factor
%
% In case the sixth input argument is provided, 
% it is interpreted as to replace the exponential in the estimate,
% i.e. it really provides an estimate of a different quantity.
% This is convenient in one of the lemmas in the Y bound.

rlam = real(lambda);
ilam = imag(lambda);
alam = abs(lambda);
rlam_f = altmid(rlam);
ilam_f = altmid(ilam);
alam_f = altmid(alam);
phi_f = altmid(phi);

if nargin<6
   expfactor=@exp;
end

Jfactor = @(rho) prefactor(rho) * expfactor(rlam_f+(sqrt(rlam_f^2*(rho+1/rho)^2+ilam_f^2*(rho-1/rho)^2))/2);

if nargin < 5 || isempty(rho0)
    % function to be optimized
    J = @(rho) Jfactor(rho) * normC0(phi_f,rho);
        
    % heuristic upper bound for rho
    rho_max = min(2 + 2*K/(alam_f), 10);

    % search for near-optimal rho
    rho0 = fminbnd(J, 1+2^-52, rho_max);

    if altisintval(phi(1))      
        rho0=intval(rho0);
    end 
end

% evaluate at near optimal rho (or at the rho given in input)
factor = prefactor(rho0) * expfactor(rlam+(sqrt(rlam^2*(rho0+1/rho0)^2+ilam^2*(rho0-1/rho0)^2))/2);
min_value = factor * normC0(phi,rho0);


