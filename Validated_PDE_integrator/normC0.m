function [ y ] = normC0(a,rho)
% normC0(a) computes an upperbound(!) for the C0-norm on [-1,1]
% hence the name of the function is a bit misleading
% Here a is assumed to be in Fourier-Chebyshev format
% normC0(a,rho) does the same but for the Bernstein-ellipse of size rho.

K = size(a,2)-1; 
k = (0:K);
if ~exist('rho','var')
    rho=1;
end

if length(rho)==1
    weights=rho.^k;
    weights=weights+1./weights;
    weights(1)=1;
    y=sum(abs(a).*weights,2);
else
    % we now set
    % weights=rho.^(0:K);
    % weights(:,2:end)=weights(:,2:end)+rho.^(-(1:K));
    % but since older Matlab versions do not support A.^B when A is a column and B a row...
    % we implement this with exp and log
    weights=exp(log(rho)*k);
    weights=weights+1./weights;
    weights(:,1)=1;
    y=sum(abs(a).*weights,2);
end
