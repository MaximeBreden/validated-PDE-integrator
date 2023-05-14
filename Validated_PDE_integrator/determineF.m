function [F,integrand,finaltime,finalboundary]=determineF(x,problem)
% determines F(x) 
% where x contains the Fourier-Chebyshev coefficients. 
% Here problem.pde contains information defining the problem
% and problem.domains contains the time steps
%
% F has the same size as x but has Fourier-chebgrid format
%
% The extra output are 
% integrand : all nonzero terms of the nonlinearity in the integrand
%             without the exponential (Fourier-Chebyshev coefficient format)
% finaltime : evaluation at the final time of integral part 
%             without initial data part and without subtracting the identity
% finalboundary : the boundary term on the final subinterval

N=(size(x,1)-1)/2;
K=size(x,2)-1;
D=size(x,3);

tau=(problem.timegrid(2:end)-problem.timegrid(1:end-1))/2; % tau=L/2

%%%%%%%%%%%%%%%%
% initial data %
%%%%%%%%%%%%%%%%

Finitial=altzeros([2*N+1,K+1,D],x(1));
p0=problem.initial.data;
% problem.initial.tailbound is not used

for d=1:D
    lambda=problem.semigroup.Lambda(:,d);
    Q=problem.semigroup.Q(:,:,d);
    Qinv=problem.semigroup.Qinv(:,:,d);
    t = cheb_grid(K,x(1));
    if d==1
        Finitial(:,:,1)= Q*(exp(tau(1)*lambda*(t+1)).*repmat(Qinv*p0,1,K+1));
    else
        pd=cheb_eval(x(:,:,d-1),1); 
        % end point of previous domain
        Finitial(:,:,d)= Q*(exp(tau(d)*lambda*(t+1)).*repmat(Qinv*pd,1,K+1));
    end
end

%%%%%%%%%%%%%%%%
% the integral %
%%%%%%%%%%%%%%%%

% the integrand (before multiplying by the exponential)
% will be fully stored (all nonzero coefficients)
Fintegral=altzeros([2*N+1,K+1,D],x(1));
maxorder=max([1,nonlinmaxdegree(problem)]);
integrand=altzeros([2*maxorder*N+1,maxorder*K+1,D],x(1));
ToepV=problem.semigroup.ToepV;

intdata=[];
for d=1:D
    lambda=problem.semigroup.Lambda(:,d);
    Q=problem.semigroup.Q(:,:,d);
    Qtau=tau(d)*Q;
    Qinv=problem.semigroup.Qinv(:,:,d);

    % extend x so that we can use convolutions to compute nonlinearities
    xx=[x(:,K+1:-1:2,d),x(:,:,d)];
                
    integranddomaind=altzeros([2*maxorder*N+1,maxorder*K+1],x(1));

    % nonlinear terms
    for p=0:length(problem.pde.polynomials)-1
        % terms with p-th derivative
        if ~isempty(problem.pde.polynomials{p+1})
            [~,phip]=convnonlinearity(xx,p,0,problem);
            phip=settensorsize(phip,[maxorder*N,maxorder*K]);
            phip=phip(:,maxorder*K+1:2*maxorder*K+1);
            phip=diag((1i.*(-maxorder*N:maxorder*N)).^p)*phip;
            integranddomaind=integranddomaind+phip;
        end
    end

    % subtract the linear semigroup term
    T=altzeros([2*N+1,2*N+1],x(1));
    DN=(1i*(-N:N));
    for p=0:length(problem.pde.polynomials)-1
        if ~isempty(problem.pde.polynomials{p+1})
            T=T+diag(DN.^p)*toeplitzshift(ToepV(:,d,p+1)); 
        end
    end
    linearpart=T*xx;
    linearpart=settensorsize(linearpart,[maxorder*N,maxorder*K]);
    % linear part has minus sign (it is subtracted since part of the semigroup)
    integranddomaind=integranddomaind - linearpart(:,maxorder*K+1:2*maxorder*K+1);

    % enforce symmetry
    if isfield(problem,'symmetry')
        if strcmp(problem.symmetry,'cosineseries')
            % cosine series: both x and integrand have real coefficients
            integranddomaind=real(integranddomaind);
        elseif strcmp(problem.symmetry,'sineseries')
            % cosine series: both x and integrand are purely imaginary
            integranddomaind=1i*imag(integranddomaind);
        end
    end

    % deal with the semigroup residue
    semigroupresidue = problem.semigroup.Residue(:,:,d)*xx;
    semigroupresidue = settensorsize(semigroupresidue,[maxorder*N,maxorder*K]);
    integranddomaind = integranddomaind - semigroupresidue(:,maxorder*K+1:2*maxorder*K+1);

    % extract the Fourier coefficients needed 
    % and multiply by Qinv 
    phi=Qinv*integranddomaind(maxorder*N+1+(-N:N),:);
           
    % compute the exponential integrals
    [intphi,intdata] = computeintegral(phi,tau(d)*lambda,K,problem,intdata);

    % multiply by Q and tau
    Fintegral(:,:,d)=Qtau*intphi;
    integrand(:,:,d)=integranddomaind;
end

F=Finitial+Fintegral;

% store the evaluation at t=1
finaltime=squeeze(Fintegral(:,1,:));
finalboundary=Finitial(:,1,D);

% subtract the identity
for d=1:D
    F(:,:,d)=F(:,:,d)-cheb_eval(x(:,:,d),t);
end

end









