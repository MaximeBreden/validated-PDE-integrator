function [ a ] = cheb_coeffs(f)
% computes the coefficients of the polynomial which
% interpolates f at the Chebyshev-points. This subroutine does not
% use the FFT (so that we can allow for an odd number interpolation
% points). 

% Initialization.
m = length(f(1,:))-1; 

if altisintval(f(1))
    ipi = intval('pi');
else
    ipi = pi;
end

% Compute 1/m * sum_j f(x_j)cos(pi*j*k/m). 
j = (0:m);
f(:,[1 m+1]) = 1/2*f(:,[1 m+1]);
a = 1/m*f*cos(ipi/m*(j'*j))'; 

% Rescale
a(:,m+1) = a(:,m+1)/2; 

end

