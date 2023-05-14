function bound = boundat1tailexp(N,d,problem)
% bounds exp(2*tau(d)*real(tlambda_n)) for all |n|>N 
% this is the exponential at t=1

tau=problem.domains/2;
P=problem.pde.order;
ToepV=problem.semigroup.ToepV;

alpha=altzeros([P/2,1],ToepV(1));

% only the even derivatives play a role
% since the odd ones lead to purely imaginary contributions to lambda
for p=0:(P/2)-1
    alpha(p+1)=(-1)^p*real(ToepV(N+1,d,2*p+1));
end
% for n>Nj the lambdas are definitely negative and decreasing
% (by simply estimating the leading term versus the sum of the others;
% this improves for the derivative)
Nalpha=sqrt(sum(abs(alpha)));
Nmax=ceil(altsup(max([Nalpha,N])));

% compute the relevant eigenvalues
nvalues=N+1:Nmax+1;
rlambda=-nvalues.^P;
if altisintval(ToepV(1))
    rlambda=intval(rlambda);
end
for p=0:(P/2)-1
    rlambda=rlambda+alpha(p+1)*nvalues.^(2*p);
end

% find the maximum
allexp = exp(2*tau(d)*rlambda);
bound=max(allexp);

end

