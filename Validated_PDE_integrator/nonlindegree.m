function l = nonlindegree(k,m,problem)
% returns the degree of the polynomial
% D_u^m g_k(u) for m=0,1 (which get k derivatives of x in the PDE)
% k will never exceed (or be equal to) problem.pde.order

if k+1>length(problem.pde.polynomials) || isempty(problem.pde.polynomials{k+1})
    % in case the term is not defined it vanishes
    l=0;
else
    l=length(problem.pde.polynomials{k+1})-m-1;
    l=max(l,0); % degree can not be negative
end
