function v=cheb_diff(u)
% determines the derivative of a Chebyshev polynomial
% The input data u is supposed to be a row vector. 
% The convention is that u = u_0 + 2*\sum u_k T_k
% The output v contains the Chebyshev coefficients of u'
    
K=size(u,2);
coeffs=2*(0:K-1)';
l=floor(K/2);
M=spdiags(repmat(coeffs,[1,l]),(1:2:K-1),K,K);
M(K,:)=[];

v=u*M';

end