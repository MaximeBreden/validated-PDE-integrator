function semigroup = fixsemigroup(x,problem)
% semigroup = fixsemigroup(y,problem)
%
% computes both the semigroup data (eigenvalues,eigenvectors,inverse)
% semigroup.Q
% semigroup.Lambda
% semigroup.Qinv
%
% of the Toeplitz/convolution matrices (one per domain)
% build up from the vectors (one per nonlinearity) 
% semigroup.ToepV
% which is real or purely imaginary if the symmetry imposes it
% 
% as well as the residue of this approximate diagonalization:
% semigroup.Residue=Q*Lambda*Qinv-Toeplitz(ToepV)
%
% These may be of interval type or float type
% depending on the type of y
% and the inverse is an outer approximation in the interval case
%
% If problem.semigroup.type='naive' then 
% Q=Qinv=identity
% Lambda=-(-N:N).^problem.pde.order
% ToepV=zero

N=(size(x,1)-1)/2;
D=size(x,3);
P=problem.pde.order;

% default values
% one for each domain
semigroup.Q=repmat(eye(2*N+1),1,1,D);
semigroup.Qinv=repmat(eye(2*N+1),1,1,D);
semigroup.Lambda=repmat(-((-N:N)').^problem.pde.order,1,D);
% and for each nonlinear term
semigroup.ToepV=zeros([2*N+1,D,P]);

% we now do the diagonalization(s)
if ~isfield(problem,'semigroup') 
    problem.semigroup.type='elaborate';
end
if ~isequal(problem.semigroup.type,'naive') 
    semigroup.type='elaborate';
    for d=1:D
        % the mean value in the domain
        average=altmid(x(:,1,d));
        for p=0:problem.pde.order-1
            % for each nonlinearity determine the linearization at the mean value
            semigroup.ToepV(:,d,p+1)=convnonlinearity(average,p,1,problem);
        end
    end
    % symmetrize
    semigroup.ToepV=...
        (semigroup.ToepV+conj(semigroup.ToepV(2*N+1:-1:1,:,:)))/2;    
    if isfield(problem,'symmetry') 
        if strcmp(problem.symmetry,'cosineseries')
            semigroup.ToepV=real(semigroup.ToepV);
        elseif strcmp(problem.symmetry,'sineseries')
            semigroup.ToepV(:,:,1:2:P)=real(semigroup.ToepV(:,:,1:2:P));
            semigroup.ToepV(:,:,2:2:P)=1i*imag(semigroup.ToepV(:,:,2:2:P));
        end
    end
    DN=(1i*(-N:N));
    for d=1:D
       % LN is the linearization at the mean value 
       LN=diag(-(abs(-N:N)').^problem.pde.order);
       for p=1:P
           LN=LN+diag(DN.^(p-1))*toeplitzshift(semigroup.ToepV(:,d,p)); 
       end
       % approximate diagonalization
       % use some perturbations to make the algorithm more stable
       perturbations=[0,1,10,100,1000]*eps;
       sizeQLQinv=zeros(size(perturbations));
       for m=1:length(perturbations)
           [mQ,mL]=eig(LN+perturbations(m));
           sizeQLQinv(m)=norm(abs(mQ)*abs(mL)*abs(inv(mQ)),1);
       end
       % pick the best one (with a penalty for the perturbation)
       [~,m0]=min(sizeQLQinv+(2*N+1)*abs(perturbations));
       [Evectors,Evalues]=eig(LN+perturbations(m0));
       semigroup.Q(:,:,d)=Evectors;
       semigroup.Lambda(:,d)=diag(Evalues);
       semigroup.Qinv(:,:,d)=inv(Evectors);
    end   
else
    semigroup.type='naive';
end

% in case of interval arithmetic, convert to intervals
if exist('intval','file') && isintval(x(1))
    semigroup.ToepV=intval(semigroup.ToepV);
    semigroup.Lambda=intval(semigroup.Lambda);
    semigroup.Qinv=intval(semigroup.Qinv);
    semigroup.Q=altzeros([2*N+1,2*N+1,D],x(1));
    for d=1:D
        % we need Q and Qinv to be exact inverses
        semigroup.Q(:,:,d)=inv(semigroup.Qinv(:,:,d));
    end
end

% now determine the residue of the approximate diagonalization
% with interval arithmetic if appropriate
semigroup.Residue=altzeros([2*N+1,2*N+1,D],x(1));

DN=(1i*(-N:N));
for d=1:D
    % LN is the linearization at the mean value 
    LN=diag(-(abs(-N:N)').^problem.pde.order);
    for p=1:P
        LN=LN+diag(DN.^(p-1))*toeplitzshift(semigroup.ToepV(:,d,p)); 
    end
    semigroup.Residue(:,:,d)=...
            semigroup.Q(:,:,d)*diag(semigroup.Lambda(:,d))*semigroup.Qinv(:,:,d)-LN;
end

end

