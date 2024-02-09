function [W,finaltime] = getWbounds(x,A,problem)
% computes the W bounds
% finaltime contains data for timestepping

%%%%%%%%%%%%%%%%%%%
%  Initialization %
%%%%%%%%%%%%%%%%%%%

N = (size(x,1)-1)/2;
K = size(x,2)-1;
D = size(x,3);

tau=problem.domains/2;
t=cheb_grid(K,x(1));
nu=problem.proof.nu;
weights=nu.^abs(-N:N)';

eps0=problem.proof.eps0;
eps1=problem.proof.eps1;
rstar=problem.proof.rstar;
vartheta=max([ones(D,1),eps0,eps1],[],2);

%%%%%%%%
% W_NK %
%%%%%%%%

WNK=altzeros([D,D,D],x(1));
P=length(problem.pde.polynomials);

normDDg=altzeros([P,D],x(1));
QDNQinvDp=altzeros([2*N+1,2*N+1,K+1,P,D],x(1));
for d=1:D
    lambda=problem.semigroup.Lambda(:,d);
    rlambda=real(lambda);
    absQ=abs(problem.semigroup.Q(:,:,d));
    absQinv=abs(problem.semigroup.Qinv(:,:,d));

    % bound on the norm in the a priori ball
    normxrstar=normL1(normC0(x(:,:,d)),nu)+vartheta(d)*rstar(d);

    Dtk=expm1div(tau(d)*rlambda*(1+t))*diag(tau(d)*(1+t));
    for p=1:P 
       % bound on norm of second derivative of nonlinearities
       normDDg(p,d)=nonlinpoly(normxrstar,p-1,2,problem,'absolute');        
       for k=1:K+1
           % construction of Xi
           QDNQinvDp(:,:,k,p,d)=absQ*diag(Dtk(:,k))*absQinv*diag(abs(-N:N).^(p-1));
       end
    end
end
        
absA=reshape(abs(A),[(2*N+1)*(K+1),D,(2*N+1)*(K+1),D]);
for d=1:D
    % permutation to allow matrix multiplication below
    Xi=permute(QDNQinvDp(:,:,:,:,d),[1 3 2 4]);
    for dd=1:D
        absAddd=reshape(absA(:,dd,:,d),(2*N+1)*(K+1),(2*N+1)*(K+1));
        normB=altzeros([P,1],x(1));
        for p=1:P
            % definition of B
            EAXi=absAddd*reshape(Xi(:,:,:,p),(2*N+1)*(K+1),2*N+1);
            % permuation to allow taking norms below
            EAXi=reshape(permute(EAXi,[1 3 2]),2*N+1,2*N+1,K+1);
            normBp=altzeros([1,K+1],x(1));
            for k=1:K+1
               normBp(k)=matrixnormL1(EAXi(:,:,k),nu);
            end
            normB(p)=normC0(normBp); 
        end
        % sum over p 
        WNK(dd,d,d)=vartheta(d)^2*sum(normB.*normDDg(:,d));
    end
end

% bounding the evalution at t=1
% factor 1/2 is taken care of in the time stepping
normQD1NQinvD=altzeros([2*N+1,P],x(1));

for p=1:P
    % evaluating Dtk at t=1
    QDat1NQinvDp=QDNQinvDp(:,:,1,p,D);
    % taking a norm
    % in the paper this operation is denoted by Upsilon
    normQD1NQinvD(:,p)=max(QDat1NQinvDp./weights,[],2);
end
% sum over p written as a matrix-matrix product
finaltime.data=normQD1NQinvD*normDDg;

%%%%%%%%%%
% WNtail %
%%%%%%%%%%

WNtail1=altzeros([D,D,D],x(1));
WNtail2=altzeros([D,D,D],x(1));
interpolK=interpolationconstant(K,0,x(1));   

% term coming from the interpolation error

for d=1:D
    lambda=problem.semigroup.Lambda(:,d);
    rlambda=real(lambda);
    absQ=abs(problem.semigroup.Q(:,:,d));
    absQinv=abs(problem.semigroup.Qinv(:,:,d));
    intlambda=expm1div(tau(d)*2*rlambda)*2*tau(d);
    D1=diag(intlambda);

    % construction of the bound
    IQLDQ=eye(2*N+1)+absQ*diag(abs(lambda))*D1*absQinv;
    normIQLDQD=altzeros([P,1],x(1));
    for p=1:P
        IQLDQDp=IQLDQ*diag(abs(-N:N).^(p-1));
        normIQLDQD(p)=matrixnormL1(IQLDQDp,nu);
    end
    % sum over p
    IQLDQDpnormDDg=sum(normIQLDQD.*normDDg(:,d));
    
    % for the second order derivative we don't see domain interaction
    WNtail1(d,d,d)=vartheta(d)^2*interpolK*tau(d)*IQLDQDpnormDDg;
end

% term coming from the improved preconditioner

for d=2:D

    QEQ = problem.interpolationerrorQEQ(:,:,d);
    
    for dd=1:D 
        Xi=permute(QDNQinvDp(:,:,:,:,dd),[1 3 2 4]);
        Areshaped=reshape(A,[(2*N+1)*(K+1),D,(2*N+1)*(K+1),D]);
        % we estimate \bar{w}^{(m-1)}(1) where m(latex)=d(code)
        % i.e. in previous domain
        Adminus1dd=reshape(Areshaped(:,d-1,:,dd),[2*N+1,K+1,2*N+1,K+1]);
        % evaluation at t=1 corresponds to sum of Chebyshev coefficients
        Eat1absAdm1dd=abs(2*sum(Adminus1dd,2)-Adminus1dd(:,1,:,:));
        Eat1absAdm1dd=reshape(Eat1absAdm1dd,2*N+1,(2*N+1)*(K+1));
        normB=altzeros([2*N+1,P],x(1));
        for p=1:P
            % take the product    
            EAXi=Eat1absAdm1dd*reshape(Xi(:,:,:,p),(2*N+1)*(K+1),2*N+1);
            % and then a norm; this is called Upsilon in the paper
            normB(:,p)=max(EAXi./weights,[],2);
        end
        % sum over p written as a matrix-vector product
        wbar=vartheta(dd)^2*normB*normDDg(:,dd);
        % we see domain interaction through the preconditioner 
        % but both derivatives need to hit the same domain
        WNtail2(d,dd,dd)=normL1(QEQ*wbar,nu);
    end
end

WNtail=WNtail1+WNtail2;

%%%%%%%%%%%%%%
% W_infinity %
%%%%%%%%%%%%%%

Winner=altzeros([D,D,D],x(1));
for d=1:D
    chi=altzeros([P,1],x(1));
    for p=1:P
        chi(p)=boundchi(p-1,N,d,problem);
    end
    Winner(d,d,d)=vartheta(d)^2*sum(chi.*normDDg(:,d));
end

% evaluation at t=1 is controled by uniform bound
finaltime.tailbound=altzeros([1,D],x(1));
for d=1:D
    finaltime.tailbound(d)=Winner(d,d,d);
end

% apply preconditioner
Ginv=estimateGinv(x,problem);
Winfinity=reshape(Ginv*reshape(Winner,D,D^2),D,D,D);

% display nothing

% final recombination of bounds
WNtail=reshape(diag(1./eps0)*reshape(WNtail,D,D^2),D,D,D);
Winfinity=reshape(diag(1./eps1)*reshape(Winfinity,D,D^2),D,D,D);
W=WNK+WNtail+Winfinity;

end

