function [Y,finaltime] = getYbounds(x,A,problem)
% computes the Y bounds
% finaltime contains data for timestepping

%%%%%%%%%%%%%%%%%%%
%  Initialization %
%%%%%%%%%%%%%%%%%%%

N = (size(x,1)-1)/2;
K = size(x,2)-1;
D = size(x,3);
dim = (2*N+1)*(K+1)*D;

%%%%%%%%%%%%%%%%
% Compute Y_KN %
%%%%%%%%%%%%%%%%

[F,integrand,finaldata,finalboundary]=determineF(x,problem);
finaltime.data=finaldata;
finaltime.boundary=finalboundary;
F=reshape(F,dim,1);
residue=reshape(A*F,[2*N+1,K+1,D]);
nu = problem.proof.nu; 
YNK=altzeros([D,1],x(1));
for d=1:D
    YNK(d) = normL1(normC0(residue(:,:,d)),nu);
end

%%%%%%%%%%%%%%%%%%
% Compute Y_infN %
%%%%%%%%%%%%%%%%%%

Ni=(size(integrand,1)-1)/2;
tau=problem.domains/2;

YNtail=altzeros([D,1],x(1));
K0 = problem.proof.interpolation.K0;
interpolK0=interpolationconstant(K0,0:K0,x(1));
prefactor = @(rho) 4*rho.^(-K0)./(rho-1);
prefactor_int = @(rho) 2*rho.^(-K0).*(rho+rho^(-1)+2)/(rho-1);
intdataK=[];
intdataK0=[];
for d=1:D
    lambda=problem.semigroup.Lambda(:,d);
    Q=problem.semigroup.Q(:,:,d);
    Qinv=problem.semigroup.Qinv(:,:,d);  
    rexplambda=expm1div(2*tau(d)*real(lambda))*2*tau(d);

    % the function appearing in the integral (apart from the exponential)
    phiY=Qinv*integrand(Ni+1+(-N:N),:,d);
    
    % first part of the interpolation bound :
    % comparing two interpolation polynomials
    t = cheb_grid(K);
    t1 = cheb_grid(K0);
    if d == 1
        tildep = Qinv*problem.initial.data;
    else
        tildep = Qinv*cheb_eval(x(:,:,d-1),1);
    end
    
    % the part corresponding to the exponentials 
    % and initial/boundary condition
    % where we postpone multiplication by tildep to reduce wrapping
    % in case the initial data are intervals (e.g. in time stepping)
    interpoly_exp = (exp(tau(d)*lambda*(t+1)));
    interpoly1_exp = (exp(tau(d)*lambda*(t1+1)));
        
    % the part corresponding to the integral
    [interpoly_int,intdataK] = computeintegral(phiY,tau(d)*lambda,K,problem,intdataK);
    [interpoly1_int,intdataK0] = computeintegral(phiY,tau(d)*lambda,K0,problem,intdataK0);
    
    % full interpolation polynomials
    % extract coefficient
    interpoly_exp = cheb_coeffs(interpoly_exp);
    interpoly_int = cheb_coeffs(interpoly_int);
    interpoly1_exp = cheb_coeffs(interpoly1_exp);
    interpoly1_int = cheb_coeffs(interpoly1_int);
    diffinterpoly_exp = interpoly1_exp - [interpoly_exp, altzeros([2*N+1,K0-K],x(1))];
    diffinterpoly_int = interpoly1_int - [interpoly_int, altzeros([2*N+1,K0-K],x(1))];
    % estimate difference
    diffinterpoly_tot = diffinterpoly_exp.*tildep + tau(d)*diffinterpoly_int;
    interperror1 = normC0(diffinterpoly_tot);   
    
    % second part of the interpolation bound :
    % interpolation errror for the higher degree polynomial
    % computed for exponential and integral terms separately

    interperror2a = altzeros([2*N+1,1],x(1));
    interperror2b = altzeros([2*N+1,1],x(1));
    for n = 1:2*N+1
        % analytic interpolation estimate for the exponential
        interperror2a(n) = analyticerrorestimate(prefactor,tau(d)*lambda(n),1,K0);
        % analytic interpolation estimate for the integral
        interperror2b(n) = tau(d)*analyticerrorestimate(prefactor_int,...
                                    tau(d)*lambda(n),phiY(n,:),K0,[],@expm1div);
    end   
    % the exponential times the inital data
    interperror2a=interperror2a.*abs(tildep); 
 
    % alternative interpolation estimate for the integral :
    % differential interpolation estimate, optimize over differentiation order q
    storeestimates = altzeros([2*N+1,K0+1],x(1)); 
    for q=0:K0
        tildeC=((tau(d)*abs(lambda)).^(q+1).*rexplambda+tau(d)*abs(tau(d)*lambda).^(q)).*normC0(phiY);
        diffcY=phiY;
        for i=1:min(q,size(phiY,2)-1) % for i larger than the degree of cY there is nothing to compute
            diffcY=cheb_diff(diffcY);
            tildeC=tildeC+tau(d)*abs(tau(d)*lambda).^(q-i).*normC0(diffcY);
        end
        storeestimates(:,q+1) = interpolK0(q+1)*tildeC;
    end
    % take the best one
    interperror2b_bis = min(storeestimates,[],2);
    interperror2b = min(interperror2b,interperror2b_bis);

    % add the two terms
    interperror2 = interperror2a+interperror2b;
    
    % full bound
    YNtail(d)=normL1(abs(Q)*(interperror1+interperror2),nu);
    
    % finally the term coming from the improved preconditioner
    if d>1
        % coupling term
        previousdomain=residue(:,:,d-1);
        wbar=cheb_eval(previousdomain,1);
        QEQ = problem.interpolationerrorQEQ(:,:,d);
        % premultiply by QEQ and take norm
        YNtail3=normL1(QEQ*abs(wbar),nu);
    else
        YNtail3=0;
    end
    
%     normL1(interperror1,nu)
%     normL1(interperror2a,nu)
%     normL1(interperror2b,nu)
%     normL1(YNtail3,nu)
% 
    YNtail(d)=YNtail(d)+YNtail3;
end

%%%%%%%%%%%%%%%%%
% Compute Y_inf %
%%%%%%%%%%%%%%%%%

Yinner=altzeros([D,1],x(1));
Fuat1=altzeros([2*Ni+1,D],x(1));
normC0F=altzeros([2*Ni+1,D],x(1));
reallambda=altzeros([2*Ni+1,D],x(1));

% integral term
for d=1:D
    lambda=-(abs(-Ni:Ni)').^problem.pde.order;
    DNi=1i*(-Ni:Ni)';
    for p=0:length(problem.pde.polynomials)-1
        if ~isempty(problem.pde.polynomials{p+1})
            lambda=lambda+(DNi).^(p)*problem.semigroup.ToepV(N+1,d,p+1);
        end
    end
    reallambda(:,d)=real(lambda);

    expterm=expm1div(2*tau(d)*reallambda(:,d))*2*tau(d); 
    % full nonlinearity
    phi=integrand(:,:,d); 
    % we want the tail only
    phi(Ni+1+(-N:N),:)=0; 
    normC0F(:,d)=normC0(phi).*expterm;

    % at the final time
    [integralat1,intdataK]=computeintegral(phi,tau(d)*lambda,1,problem,intdataK);
    Fuat1(:,d)=tau(d)*integralat1(:,1);
end

for d=1:D
    Yinner(d)=normL1(normC0F(:,d),nu);
end

% tailbound at t=1
finaltime.tailbound=normL1(Fuat1,nu);

% C0 estimate decay of the tailbound (initial data term)
normC0Finitial=boundC0tailexp(N,1,problem)*problem.initial.tailbound;

% estimate decay of the tailbound (initial data term) at t=1
Fat1initial=boundat1tailexp(N,1,problem)*problem.initial.tailbound;

% apply preconditioner
Ginv=estimateGinv(x,problem);

% Missing from Yinner(1) is the initial value term at t=1.
% It does not estimate Fat1initial by normC0Finitial
% although the integral part is estimated uniformly.
% Of course it only makes a difference in the first domain
% and only if problem.initial.tailbound is nonzero
ivpcorrection=[Fat1initial;altzeros([D-1,1],x(1))];

Yinfinity=Ginv*(Yinner+ivpcorrection);

% Override the first domain
% where we need to add the C0 norm of the initial value term
Yinfinity(1)=Yinner(1)+normC0Finitial;

% display
disp(['size of YNK is ',num2str(altsup(max(YNK)))]);
disp(['size of YNtail is ',num2str(altsup(max(YNtail)))]);
disp(['size of Yinfinity is ',num2str(altsup(max(Yinfinity)))]);

% final recombination of bounds
eps0 = problem.proof.eps0;
eps1 = problem.proof.eps1;
Y = altsup(YNK+YNtail./eps0+Yinfinity./eps1);

end

