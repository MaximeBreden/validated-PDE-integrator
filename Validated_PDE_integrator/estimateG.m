function interpolationerrorQEQ = estimateG(x,problem)
% computes the interpolation error on the exponential functions 
% needed to bound G

N = (size(x,1)-1)/2;
K = size(x,2)-1;
D = size(x,3);
tau = problem.domains/2;

K0 = problem.proof.interpolation.K0;
prefactor = @(rho) 4*rho.^(-K0)./(rho-1);
interpolation_error = altzeros([2*N+1,1],x(1));
interpolationerrorQEQ = altzeros([2*N+1,2*N+1,D],x(1));

for d = 1:D
    lambda = problem.semigroup.Lambda(:,d);
    
    % first part of the bound:
    % the difference between two interpolation polynomials
    t = cheb_grid(K);
    t1 = cheb_grid(K0);    
    InterpExp = cheb_coeffs(exp(tau(d)*lambda*(t+1)));
    InterpExp1 = cheb_coeffs(exp(tau(d)*lambda*(t1+1)));
    diff_poly_exp = InterpExp1 - [InterpExp, altzeros([2*N+1,K0-K],x(1))];

    error_diff_poly_exp = normC0(diff_poly_exp);
    
    % second part of the bound: 
    % interpolation errror for the higher degree polynomial
    for n = 1:2*N+1
        interpolation_error(n) = analyticerrorestimate(prefactor,tau(d)*lambda(n),1,K0);
    end
    interpolation_error = interpolation_error + error_diff_poly_exp;
    
    % multiply by Q and Qinv
    Q=problem.semigroup.Q(:,:,d);
    Qinv=problem.semigroup.Qinv(:,:,d);
    interpolationerrorQEQ(:,:,d) = abs(Q)*diag(interpolation_error)*abs(Qinv);
end


end