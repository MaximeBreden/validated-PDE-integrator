function t = cheb_grid(K,intvaltest)
% returns the Chebyshev grid with K+1 grid points
% if intvaltest is true t is interval-valued

if exist('intvaltest','var') && altisintval(intvaltest)
   ipi = intval('pi'); 
else
   ipi = pi; 
end

t = cos((0:K)*ipi/K);
    
end

