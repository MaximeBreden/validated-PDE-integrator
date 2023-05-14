function [ y ] = cheb_eval(a,t)
% evaluates the Chebyshev polynomial 
% corresponding to the Chebyshev coefficients a

K=size(a,2)-1;
k=(0:K)';
cosktheta=cos(k*acos(t));

a(:,2:end)=2*a(:,2:end);
% absorb the factor 2 in the definition of the Chebyshev coefficients 
y=a*cosktheta;

end
