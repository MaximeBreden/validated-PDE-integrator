function weights = cheb_quadrature_weights(K,intvaltest)

% Computes the Clenshaw-Curtis quadrature weights, by applying the DCT to
% the vector b containing the integrals of the Chebyshev polynomials over
% [-1,1].
% If intvaltest is true then the weights are interval-valued

if exist('intvaltest','var') && altisintval(intvaltest)
    ind = intval(0:K);
    ipi = intval('pi');
else
    ind = 0:K;
    ipi = pi;
end

% the Chebyshev nodes are given by cos(theta)
theta = (K-ind')*ipi/K; 

% M is the transpose of the matrix sending a vector of values at Chebyshev 
% nodes to a vector of Chebyshev coefficients. That is:
% if u = [f(x_0);...;f(x_K)], and v = [f_0;...;f_K], where the f_k are such
% that f = f_0*T_0 + ... + f_K*T_K, the T_k being the Chebyshev polynomials
% of the first kind, then transpose(M)*u = v.
M =  cos( theta * ind) / K; 
M([1,K+1],:) = M([1,K+1],:) / 2;
M(:,2:K,:) = 2*M(:,2:K); 

% b = [b_0;...;b_K], where b_k = int_{-1}^1 T_k(x)dx
b = 2./(1-ind'.^2);
b(2:2:end) = 0;

weights = M*b;