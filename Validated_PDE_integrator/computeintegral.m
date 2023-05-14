function [I,integrationdata,errordata] = ...
                    computeintegral(phi,lambda,KK,problem,integrationdata,errordata)

% This subroutine computes the integrals
% 
%   int_{-1}^{t} exp(lambda*(t-s))*phi(s) ds 
%
% for the Chebyshev nodes t=t(k), k=1..KK 
% where phi is given in terms of Chebyshev coefficients
% (size 2*Nphi+1 by Kphi+1). 
%
% The integrationdata contains precomputed values of various constants.
% It gets updated depending on phi and lambda and KK
% If integrationdata==[] then all data is computed
%
% The last argument errordata is optional. 
% If given, it forces the error bound to use 
% the Bernstein ellipse of size errordata.rho 
% (unless lambda has been changed)
% If not given, an optimal rho is computed inside analyticerrorestimate.
%
% Note that KK and Kphi may be different from K, Nphi may be different from N
% The length of lambda is assumed to be 2*Nphi+1
% and if phi is intval then lambda is assumed to be intval as well

%%%%%%%%%%%%%%%%%%
% Initialization %
%%%%%%%%%%%%%%%%%%

if (isfield(problem,'type') && strcmp(problem.type,'proof')) || ...
                                                  altisintval(phi(1))
    % proof
    K1=problem.proof.integrals.K1;
    errorestimate=true;
    if altisintval(phi(1))
        % with interval arithmetic
        useintervals = true;
        fraction = intval(64)/intval(15);
    else
        useintervals = false;
        fraction=64/15;
    end
else
    % numerics only
    K1=problem.numerics.integrals.K1;
    errorestimate=false;
    useintervals=false;
end

if isempty(integrationdata)
    updateall=true;
elseif ~(altisintval(integrationdata.grids{1}(1))==useintervals)
    % need to change from floats to intvals or vice versa
    updateall=true;
else
    updateall=false;
    if ~(integrationdata.cheb.K1==K1)
        updateK1=true;
    else
        updateK1=false;
    end
    if ~(length(integrationdata.grids)==KK+1)
        updateKK=true;
    else
        updateKK=false;
    end
    if ~(length(integrationdata.lambda)==length(lambda)) || ...
                                    ~all(integrationdata.lambda==lambda)
        updatelambda=true;
    else
        updatelambda=false;
    end
    if ~(size(integrationdata.phigridmatrix{1},1)==size(phi,2))
        updateKphi=true;
    else
        updateKphi=false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update what needs to be updated %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kphi=size(phi,2)-1;
Nphi=(size(phi,1)-1)/2;

if updateall || updateK1
    integrationdata.cheb=fixcheb(K1,phi(1));
end
KK1=integrationdata.cheb.KK1;
wK1=integrationdata.cheb.weights;
tK1=integrationdata.cheb.grid;
t = cheb_grid(KK,phi(1));

if updateall || updateK1 || updateKK    
    integrationdata.grids=cell(KK+1,1);
    for k=1:KK+1
        % rescaled quadrature grid for the exponential (which has factor t-s)
        integrationdata.grids{k}=(t(k)+1)/2*(1-tK1);
    end
end
grids=integrationdata.grids;

if updateall || updateK1 || updateKK || updatelambda
    integrationdata.explambda=cell(KK+1,1);
    for k=1:KK+1
        % evaluate the exponential on quadrature grid
        integrationdata.explambda{k}=exp(lambda*grids{k});
    end
    integrationdata.lambda=lambda;
    lambdasupdated=true;
else
    lambdasupdated=false;
end
explambda=integrationdata.explambda;

if updateall || updateK1 || updateKK || updateKphi
    integrationdata.phigridmatrix=cell(KK+1,1);
    for k=1:KK+1
        % rescaled quadrature grid for the integrand without exponential
        grid_phi = (t(k)+1)/2*(tK1+1)-1; 
        % matrix corresponding to evaluating Chebyshev polynomials on quadrature grid
        integrationdata.phigridmatrix{k}=cheb_evalmatrix(grid_phi,Kphi);
    end
end
phigridmatrix=integrationdata.phigridmatrix;

%%%%%%%%%%%%%%
% Quadrature %
%%%%%%%%%%%%%%

integrand=altzeros([2*Nphi+1,KK+1,KK1+1],phi(1));
for k=1:KK+1
    % evaluate the Chebyshev polynomial phi and the integand on quadrature grid
    phiongrid = phi*phigridmatrix{k};% 
    integrand(:,k,:) = explambda{k}.*phiongrid;
end
% multiply the integrand at the quadrature grid by the quadrature weights
integrand=reshape(integrand,(2*Nphi+1)*(KK+1),KK1+1);
quadrature = integrand*wK1;
quadrature = reshape(quadrature,2*Nphi+1,KK+1);

I = quadrature*diag((t+1)/2);

if ~errorestimate
    errordata=[];
    return
end

%%%%%%%%%%%%%%%%%%
% Error estimate %
%%%%%%%%%%%%%%%%%%

if ~exist('errordata','var') || lambdasupdated
    optimizerho=true;
else
    optimizerho=false;
    errorfactor=errordata.factor;
    rho=errordata.rho;
end

% bound the quadrature error
r = altzeros([2*Nphi+1,KK+1],phi(1));
if optimizerho
    prefactor = @(x) x^(-KK1)/(x^2-1); % KK1 is even
    errordata.rho=altzeros([2*Nphi+1,1],phi(1));
    errordata.factor=altzeros([2*Nphi+1,1],phi(1));
    for n=1:length(lambda)
        [err,rho0,factor] =  analyticerrorestimate(prefactor,lambda(n),phi(n,:),KK1);
        % include the factor 64/15 and (t+1)/2 missing from the interpolation estimate
        r(n,:) = fraction*err*(t+1)/2;
        errordata.rho(n)=rho0;
        errordata.factor(n)=factor;
    end
    % update 
    errordata.lambda=lambda;
else
    err = errorfactor .* normC0(phi,rho);
    r = fraction*err*(t+1)/2;
end

if useintervals
    I=I+ cintval(0,altsup(r));
else 
    %non-rigorous, no intervals
    I=I+sign(I).*r;
end

% bound of the whole integral 
% for large eigenvalues this might be better
% than quadrature + quadrature error 
upperbound = normC0(phi)*(1+t).*expm1div(real(lambda)*(1+t));
if useintervals
    I=intersect(I,cintval(0,altsup(upperbound))); 
else 
    %non-rigorous, no intervals
    I=sign(I).*min(abs(I),upperbound);
end

end




