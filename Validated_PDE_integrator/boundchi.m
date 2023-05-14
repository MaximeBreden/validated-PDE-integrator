function chi = boundchi(j,N,d,problem)
% computes the bound chi needed to estimate Fourier truncation errors
%
% bounds n^j (exp(2*tau(d)*real(tlambda_n))-1)/real(tlambda_n) for all |n|>N

tau=problem.domains/2;
P=problem.pde.order;
ToepV=problem.semigroup.ToepV;

alpha=altzeros([P/2,1],ToepV(1));
djalpha=altzeros([P/2,1],ToepV(1));

% only the even derivatives play a role
% since the odd ones lead to purely imaginary contributions to lambda
for p=0:(P/2)-1
    alpha(p+1)=(-1)^p*real(ToepV(N+1,d,2*p+1));
    djalpha(p+1)=abs(alpha(p+1))*abs(2*p-j)/(P-j);
end
Ndjalpha=sqrt(sum(djalpha));
% beyond Nmax, n^j/real(tlambda_n) is an upper bound and is decreasing
Nmax=ceil(altsup(max([Ndjalpha,N+1])));

% compute the relevant eigenvalues
nvalues=N+1:Nmax;
rlambda=-nvalues.^problem.pde.order;
for p=0:(problem.pde.order/2)-1
    rlambda=rlambda+alpha(p+1)*nvalues.^(2*p);
end

% find the maximum
njexp = nvalues.^j.*expm1div(2*tau(d)*rlambda)*2*tau(d);
tildechi=max(njexp);

% check whether the n^j/real(tlambda_n) upper bound is already good enough,
% or whether we need to compute a bit further
tNmax = Nmax;
rlambda = rlambda(end);
while tNmax^j/(-rlambda) > tildechi
    tNmax = tNmax+1;
    rlambda = -tNmax^problem.pde.order;
    for p=0:(problem.pde.order/2)-1
        rlambda=rlambda+alpha(p+1)*tNmax^(2*p);
    end
end
if tNmax == Nmax % the n^j/real(tlambda_n) upper bound was already good enough
    chi = tildechi;
else
    % compute the missing eigenvalues
    nvalues=Nmax+1:tNmax;
    rlambda=-nvalues.^problem.pde.order;
    for p=0:(problem.pde.order/2)-1
        rlambda=rlambda+alpha(p+1)*nvalues.^(2*p);
    end

    % find the maximum
    njexp = nvalues.^j.*expm1div(2*tau(d)*rlambda)*2*tau(d);
    chi = max(tildechi,max(njexp));
end

end
