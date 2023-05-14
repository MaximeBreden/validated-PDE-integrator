function [DF,Dphifull]=determineDF(x,problem)
% determines DF(x) 
% where x contains the Fourier-Chebyshev coefficients. 
% Here problem.pde contains information defining the problem
% and problem.domains contains the time steps
%
% DF is a tensor of size (2N+1)x(K+1)xDx(2N+1)x(K+1)xD
%
% The extra output Dphifull is contains all nonzero terms of the derivatives 
% of the nonlinearities in the integrand without the exponential 
% and without the semigroup residue
% (Fourier-Chebyshev coefficient format)

N=(size(x,1)-1)/2;
K=size(x,2)-1;
D=size(x,3);

tau=(problem.timegrid(2:end)-problem.timegrid(1:end-1))/2; % tau=L/2

%%%%%%%%%%%%%%%%%%
% initialization %
%%%%%%%%%%%%%%%%%%

t = cheb_grid(K);
DF=altzeros([2*N+1,K+1,D,2*N+1,K+1,D],x(1));
% matrix used for the derivative of the identity
% which we interpret to go from Cheb coefficients to values at Cheb nodes
% transposed because of multiplication on the right
evalmatrix=cheb_evalmatrix(t,K)';
minusidentity=-reshape(evalmatrix,[1,K+1,1,1,K+1,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparation for the integral %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxsize=max([1,nonlinmaxdegree(problem)])+1;
KK=(maxsize-1)*K;   

Dphifull=cell(length(problem.pde.polynomials),1);
for p=0:length(problem.pde.polynomials)-1
    degree=max([1,nonlindegree(p,1,problem)]);
    Dphifull{p+1} = altzeros([2*degree*N+1,2*degree*K+1,D],x(1));
end
    
ToepV=problem.semigroup.ToepV;

intdata=[];
errordata=[];
for d=1:D
    lambda=problem.semigroup.Lambda(:,d);
    Q=problem.semigroup.Q(:,:,d);
    Qtau=tau(d)*Q;
    Qinv=problem.semigroup.Qinv(:,:,d);
    
    % extend x so that we can use convolutions to compute nonlinearities
    xx=[x(:,K+1:-1:2,d),x(:,:,d)];

    % linear terms
    T=altzeros([2*N+1,2*N+1],x(1));
    DN=(1i*(-N:N));
    for p=0:length(problem.pde.polynomials)-1  
        if ~isempty(problem.pde.polynomials{p+1})
            T=T+diag(DN.^(p))*toeplitzshift(ToepV(:,d,p+1)); 
        end
    end 
    semigroupresidue = problem.semigroup.Residue(:,:,d);

    % nonlinear terms 
    Dphi=cell(length(problem.pde.polynomials),1);
    for p=0:length(problem.pde.polynomials)-1
        if ~isempty(problem.pde.polynomials{p+1})
            % terms with p-th derivative
            [~,Dphip]=convnonlinearity(xx,p,1,problem);
            Dphifull{p+1}(:,:,d)=Dphip;
            Dphi{p+1}=settensorsize(Dphip,[2*N,maxsize*K]);
        end
    end

    reoptimizeerrorbound=true;
    for n=-N:N
        for l=0:K
            % derivative w.r.t n-the Fourier and l-th Chebyshev mode
            % by shift of the appropriate Fourier-Cheybshev block

            Dphiall=altzeros([2*N+1,2*KK+1],x(1));
            for p=0:length(problem.pde.polynomials)-1
                if ~isempty(problem.pde.polynomials{p+1})
                    Dp=Dphi{p+1}(2*N+1+(-N:N)-n,maxsize*K+1+(-KK:KK)-l);
                    if l>0
                        % Chebyshev is a cosine series, i.e. add shift over -l
                        Dp=Dp+Dphi{p+1}(2*N+1+(-N:N)-n,maxsize*K+1+(-KK:KK)+l);
                    end
                    Dphiall=Dphiall+diag((1i*(-N:N)).^p)*Dp;
                end
            end
            Dphiall=Dphiall(:,KK+1:2*KK+1);

            % add derivative of (minus) linear term
            Dphiall(:,l+1)=Dphiall(:,l+1)-T(:,N+1+n);
            if isfield(problem,'symmetry') && contains(problem.symmetry,'sineseries') 
                % enforce symmetry in case of cosine or sine series
                Dphiall=real(Dphiall);
            end
            % deal with the semigroup residues
            Dphiall(:,l+1)=Dphiall(:,l+1)-semigroupresidue(:,N+1+n);

            Dphiall=Qinv*Dphiall;

            % compute the exponential integrals
            if reoptimizeerrorbound
                % We only go here for the first of the derivatives 
                % i.e. if l==0 && n==-N
                % To get optimal bounds one should optimize for all l and n
                % but that is slow
                [intDphi2,intdata,errordata] = ...
                    computeintegral(Dphiall,tau(d)*lambda,K,problem,intdata);
                reoptimizeerrorbound=false;
            else
                [intDphi2,intdata] = ...
                  computeintegral(Dphiall,tau(d)*lambda,K,problem,intdata,errordata);
            end
            DF(:,:,d,n+N+1,l+1,d)=Qtau*intDphi2;   
        end % end of loop over Chebyshev mode derivative
    end % end of loop over Fourier mode derivative
    
    for n=1:2*N+1
        % add the derivative of minus identity
        DF(n,:,d,n,:,d)=DF(n,:,d,n,:,d) + minusidentity;
    end
    if d>1
        % derivative of the boundary term
        % this is the only way different domains are coupled
        for k=1:K+1    
            % since u(1) = a_0+ 2*sum_{k>1} a_k in terms of Chebychev coefficients
            rightend=2*repmat(Q*diag(exp(tau(d)*lambda*(t(k)+1)))*Qinv,1,1,K+1);
            rightend(:,:,1)=rightend(:,:,1)/2;
            DF(:,k,d,:,:,d-1)=reshape(rightend,[2*N+1,1,1,2*N+1,K+1,1]);
        end
    end
end % end of loop over domain

end









