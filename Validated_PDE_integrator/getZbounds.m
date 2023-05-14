function [Z,finaltime] = getZbounds(x,A,problem)
% computes the Z bounds
% finaltime contains data for timestepping

%%%%%%%%%%%%%%%%%%%
%  Initialization %
%%%%%%%%%%%%%%%%%%%

N = (size(x,1)-1)/2;
K = size(x,2)-1;
D = size(x,3);
dim = (2*N+1)*(K+1)*D;

%%%%%%%%%%%%%%%%%%%%%
% Compute ZNK bound %
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% Compute Z0 %
%%%%%%%%%%%%%%

nu=problem.proof.nu;

% compute I-A*DF
Z0=altzeros([D,D],x(1));
[DF,Dphifull]=determineDF(x,problem);
DF=reshape(DF,dim,dim);
IADFNK=abs(eye(dim)-A*DF);
IADFNK=reshape(IADFNK,[2*N+1,K+1,D,2*N+1,K+1,D]);

% compute the l1 norms 
weights=repmat(nu.^abs(-N:N)',1,K+1);
weights(:,2:K+1)=2*weights(:,2:K+1);
weights=reshape(weights,1,(2*N+1)*(K+1));
for d=1:D
    for dd=1:D
        IAFNKddd=reshape(IADFNK(:,:,d,:,:,dd),(2*N+1)*(K+1),(2*N+1)*(K+1));
        Z0(d,dd)=max((weights*abs(IAFNKddd))./weights);
    end
end

%%%%%%%%%%%%%%
% Compute Z1 %
%%%%%%%%%%%%%%

% initializations
Z1=altzeros([D,D],x(1));
DFdiff=altzeros([2*N+1,K+1,D],x(1));
t=cheb_grid(K);
eps0=problem.proof.eps0;
eps1=problem.proof.eps1;
vartheta=max([ones(D,1),eps0,eps1],[],2);
tau=problem.domains/2;
P=length(problem.pde.polynomials);

% used to store some data needed later
Phitilde=altzeros([P,D],x(1));  
finaltime.data=altzeros([2*N+1,D],x(1));

for d=1:D
    lambda=problem.semigroup.Lambda(:,d);
    rlambda=real(lambda);
    absQ=abs(problem.semigroup.Q(:,:,d));
    Qinv=problem.semigroup.Qinv(:,:,d);
    absQinv=abs(Qinv);

    % size of full convolution derivatives
    Np=zeros([P,1]);
    Kp=zeros([P,1]);
    for p=1:P
        Np(p)=(size(Dphifull{p}(:,:,d),1)-1)/2;
        Kp(p)=(size(Dphifull{p}(:,:,d),2)-1)/2;
    end
    Nmax=max(Np);
    Kmax=max(Kp);

    DjPhi=altzeros([2*N+1,1],x(1));
    DjPhicheck=altzeros([2*N+1,1],x(1));
    for p=1:P
        Vp=problem.semigroup.ToepV(:,d,p);
        % nonlinear term
        Dphip=settensorsize(Dphifull{p}(:,:,d),[Nmax,Kmax]);
        Dphip=Dphip(:,Kmax+1:2*Kmax+1);
        Dphipmod=Dphip;
        Dphipmod(Nmax+1+(-N:N),1)=Dphipmod(Nmax+1+(-N:N),1)-Vp;
        
        % get the bound on this term
        [Phi,Phicheck,Phitilde(p,d)]=...
               boundphi(normC0(Dphipmod),normC0(Dphip),N,eps0(d),eps1(d),nu);

        %collect all of them
        DjPhi=DjPhi+diag((abs(-N:N)).^(p-1))*Phi;
        DjPhicheck=DjPhicheck+diag((abs(-N:N)).^(p-1))*Phicheck;
    end

    % diagonalization residue term
    absQinvR=abs(Qinv*problem.semigroup.Residue(:,:,d));
    Ypsilon=max(absQinvR*diag(1./nu.^abs(-N:N)),[],2);

    for k=1:K+1
        intlambda=expm1div(tau(d)*(1+t(k))*rlambda)*(1+t(k));
        Dtk=diag(intlambda);
        % This is D_N(t) except for a missing factor L/2 
        % which we correct for below
        DFdiff(:,k,d)=tau(d)*absQ*Dtk*(absQinv*DjPhi+eps0(d)*Ypsilon);
        if k==1 
           finaltime.data(:,d)=tau(d)*absQ*Dtk*...
                          (absQinv*DjPhicheck+max(1,eps0(d))*Ypsilon);
        end
    end

end

% premultiply by A
DFdiff=reshape(DFdiff,(2*N+1)*(K+1),D);
absA=reshape(abs(A),[(2*N+1)*(K+1),D,(2*N+1)*(K+1),D]);
for d=1:D
    DFdiffd=DFdiff(:,d);
    for dd=1:D
       absAdd=reshape(absA(:,dd,:,d),(2*N+1)*(K+1),(2*N+1)*(K+1));
       absAddDFdiffd=reshape(absAdd*DFdiffd,2*N+1,K+1);
       Z1(dd,d)=normL1(normC0(absAddDFdiffd),nu);
    end
end

ZNK=Z0+Z1;

%%%%%%%%%%
% ZNtail %
%%%%%%%%%%

ZNtail1=altzeros([D,D],x(1));
ZNtail2=altzeros([D,D],x(1));
ZNtail3=altzeros([D,D],x(1));
interpolK=interpolationconstant(K,0,x(1));    

%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolation estimate %
%%%%%%%%%%%%%%%%%%%%%%%%%%

for d=1:D
    lambda=problem.semigroup.Lambda(:,d);
    rlambda=real(lambda);
    Q=problem.semigroup.Q(:,:,d);
    Qinv=problem.semigroup.Qinv(:,:,d);

    % this is D_N(1)
    intlambda=expm1div(tau(d)*2*rlambda)*2*tau(d);
    D1=diag(intlambda);
    term1=0;
    IQLDQ=eye(2*N+1)+abs(Q)*diag(abs(lambda))*D1*abs(Qinv);
    % terms corresponding to interpolation error in integral terms
    for p=1:P
        IQLDQD=IQLDQ*diag(abs(-N:N).^(p-1));
        normIQLDQD=matrixnormL1(IQLDQD,nu);
        term1=term1+normIQLDQD*Phitilde(p,d);
    end
        
    % term corresponding to diagonalization residue
    IQLDQ=eye(2*N+1)+abs(Q)*diag(abs(lambda))*D1*abs(Qinv);
    normIQLDQ=matrixnormL1(IQLDQ,nu); 
    Rnorm=matrixnormL1(absQinvR,nu);
    term2=max(1,eps0(d))*normIQLDQ*Rnorm;
    
    ZNtail1(d,d)=interpolK*tau(d)*(term1+term2);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% preconditioner terms %
%%%%%%%%%%%%%%%%%%%%%%%%

for d=2:D
    
    QEQ = problem.interpolationerrorQEQ(:,:,d);
    
    % two terms coming from the improved preconditioner
    
    for dd=1:D
       % take the matrix mapping to the prior domain
       % as that one is showing up as initial data in the current domain
       IADFNKdm1dd=reshape(IADFNK(:,:,d-1,:,:,dd),[2*N+1,K+1,2*N+1,K+1]);
       % evaluate at t=1 (sum of Chebyshev coefficients)
       at1=reshape(2*sum(IADFNKdm1dd,2)-IADFNKdm1dd(:,1,:,:),[2*N+1,(2*N+1)*(K+1)]);
       % interpret as n maps from Fourier-Cheb coefficient space to scalar
       % and compute norm
       wbar=max(abs(at1)./weights,[],2);
       % premultiply by QEQ and take norm
       ZNtail2(d,dd)=normL1(QEQ*wbar,nu);
    end
    
    for dd=1:D
       % take the matrix mapping to the previous domain (as before)
       absAdm1=reshape(absA(:,d-1,:,dd),(2*N+1)*(K+1),(2*N+1)*(K+1));
       absADFdiff=reshape(absAdm1*DFdiff(:,dd),2*N+1,K+1);
       % evaluate at t=1
       wbar=cheb_eval(absADFdiff,1);
       % premultiply by QEQ and take norm
       ZNtail3(d,dd)=normL1(QEQ*wbar,nu);
    end

end

ZNtail=ZNtail1+ZNtail2+ZNtail3;

%%%%%%%%%%%%%%
% Z_infinity %
%%%%%%%%%%%%%%

Zinner=altzeros([D,D],x(1));

for d=1:D
    rlambdaN1=-(N+1)^problem.pde.order;
    normgtilde=altzeros([P,1],x(1));

    for p=1:P
        Dphip=Dphifull{p}(:,:,d);
        Kp=(size(Dphip,2)-1)/2;
        % full derivative
        gtilde=Dphip(:,Kp+1:2*Kp+1);
        Vp=problem.semigroup.ToepV(:,d,p);
        Np=(size(Dphifull{p}(:,:,d),1)-1)/2;
        % subtract the average
        gtilde(Np+1,1)=gtilde(Np+1,1)-Vp(N+1);
        normgtilde(p)=normL1(normC0(gtilde),nu);
        rlambdaN1=rlambdaN1+(1i*(N+1))^(p-1)*Vp(N+1);
    end

    chi=altzeros([P,1],x(1)); 
    for p=1:P
        chi(p)=boundchi(p-1,N,d,problem);
    end
 
    Zinner(d,d)=vartheta(d)*sum(chi.*normgtilde);
end

finaltime.tailbound=diag(Zinner)';

% apply preconditioner
Ginv=estimateGinv(x,problem);
Zinfinity=Ginv*Zinner;

%display
disp(['size of ZNK is ',num2str(altsup(max(diag(ZNK))))]);
disp(['size of ZNtail is ',num2str(altsup(max(diag(ZNtail))))]);
disp(['size of Zinfinity is ',num2str(altsup(max(diag(Zinfinity))))]);

% final recombination of bounds
Z=ZNK+diag(1./eps0)*ZNtail+diag(1./eps1)*Zinfinity;

end
