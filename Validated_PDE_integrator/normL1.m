function [ y ] = normL1(a,nu)
% normL1(a,nu) computes the l^1_nu norm of a vector a
% the default value for nu is 1
% if a is a matrix the norm is computed columnwise

if ~exist('nu','var')
    nu=1;
end
N = (size(a,1)-1)/2; 
weights = nu.^abs(-N:N)'; 
y = sum(abs(a).*weights,1);

end
