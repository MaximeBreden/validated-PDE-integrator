 function [success,rsol,finaltime,bounds,errorbound]=provesol(x,problem)
% Takes the approximate solution x and problem parameters and 
% attempts to prove a solution closeby exists
% If x consists of intervals all computations will be rigorously verified
% using interval arithmetic in Intlab 
% (otherwise it is based on floating point arithmetic)
% Here x is an (2N+1)x(K+1)xD tensor
% and we assume x is already conjugate symmetric in the Fourier component
%
% If the radii polynomial is negative for some positive r, 
% then success=true and rsol are the succesfull radii 
%
% finaltime contains bounds on some integral terms to enclose the solution 
% at the end of the time interval for time stepping
%
% Additionally, bounds contains some additional data about the bounds 
% (for saving) from which the radii polynomials have been built 
%
% If the method is unsuccesfull then success=0 and rsol=NaN

success=false;
rsol=NaN;

N = (size(x,1)-1)/2;
K = size(x,2)-1;
D = size(x,3);
dim = (2*N+1)*(K+1)*D;

if altisintval(x(1)) && ~altisintval(problem.semigroup.ToepV(1))
    % if data of semigroup is not of interval type
    % recompute semigroup to get exact inverses Q and Qinv
    problem.semigroup = fixsemigroup(x,problem);
end
    
% determine A
disp('Defining the matrix A')

xfloat=altmid(x);
problemfloat=problem;
problemfloat.type='numerics';
% make sure these are all floats
problemfloat.semigroup.ToepV=altmid(problem.semigroup.ToepV);
problemfloat.semigroup.Q=altmid(problem.semigroup.Q);
problemfloat.semigroup.Qinv=altmid(problem.semigroup.Qinv);
problemfloat.semigroup.Lambda=altmid(problem.semigroup.Lambda);
problemfloat.semigroup.Residue=altmid(problem.semigroup.Residue);
problemfloat.timegrid=altmid(problemfloat.timegrid);

DF=determineDF(xfloat,problemfloat);
DF=reshape(DF,dim,dim);
if isfield(problem,'symmetry') && ...
        (strcmp(problem.symmetry,'cosineseries') || strcmp(problem.symmetry,'sineseries'))
    % the Jacobian is (approximately) real in case of cosine or sine series
    DF=real(DF);
end
A=inv(DF);

% symmetrize A (not required for the proof)
A0=reshape(A,[2*N+1,K+1,D,2*N+1,K+1,D]);
for d=1:D
    % since DF is lower triangular so should A
    A0(:,:,d,:,:,d+1:D)=0;
end
A1=conj(flip(flip(A0,1),4));
A=reshape((A0+A1)/2,[dim,dim]);

disp(['norm of A is ',num2str(norm(A,1))]);
sizeA=size(A,1);
disp(['dimension of A is ',int2str(sizeA)]);

% Now the interval arithmetic computations start
disp('Starting the proof')
problem.type='proof';
% turn all data into intervals if x is intval
if altisintval(x(1))
    disp('Including interval arithmetic')
    A=intval(A); 
    problem.timegrid=intval(problem.timegrid);
    problem.initial.tailbound=intval(problem.initial.tailbound);
    problem.proof.nu=intval(problem.proof.nu);
    problem.proof.eps0=intval(problem.proof.eps0);
    problem.proof.eps1=intval(problem.proof.eps1);
    problem.proof.rstar=intval(problem.proof.rstar);
	if ~isintval(problem.initial.data)
		initialdataradii=1.414214*eps(problem.initial.data);
		problem.initial.data=midrad(problem.initial.data,initialdataradii);     
	end
else
    disp('No intervals yet');
    problem.initial.data=altmid(problem.initial.data);
    problem.timegrid=altmid(problem.timegrid);
end
problem.domains=(problem.timegrid(2:end)-problem.timegrid(1:end-1));

% degree of the higher order polynomial which we use to approximate the interpolation error
problem.proof.interpolation.K0 = max(problem.proof.interpolation.K0,K);

% the bound on the preconditioner introduced for the domain decomposition
interpolationerrorQEQ = estimateG(x,problem);
problem.interpolationerrorQEQ = interpolationerrorQEQ;

% get the Y bound
[Y,finaltime.Y] = getYbounds(x,A,problem);
if isfield(problem,'displaybounds') && problem.displaybounds
    disp('The Y-bound:')
    disp(num2str(altsup(Y)));    
end
disp(['Rough norm of the Y-bound: ',num2str(altsup(max(Y)))]);

% get the Z bound
[Z,finaltime.Z] = getZbounds(x,A,problem);
if isfield(problem,'displaybounds') && problem.displaybounds
    disp('The Z-bound:');
    disp(num2str(altsup(Z)));  
end
disp(['Dominant eigenvalue of the Z-matrix is roughly ',num2str(eigs(altsup(Z),[],1))]); 

% get the W bound
[W,finaltime.W] = getWbounds(x,A,problem);
% size of W bound: sum(pagenorm(W,Inf)
sizeW=sum(max(sum(altsup(W),2),[],1,'includenan'));
disp(['Rough size of W: ',num2str(sizeW)]); 

% find a point where the radii polynomials are negative
[rminvec,eta]=polynomialsnegative(Y,Z,W);
rstarvec=problem.proof.rstar;

if any(isnan(rminvec)) || ~all(rminvec>0)
    % failure :-(
    disp('Failed due not finding a positive rmin');
elseif all(rminvec <= rstarvec)
    % success :-)
    success=true;  
    rsol=rminvec;
else
    % failure :-(
    disp('Failed due to rmin not being smaller than rstar');
    disp(['rmin  = [' ,num2str(rminvec'),' ]']); 
    disp(['rstar = [' ,num2str(altsup(rstarvec')),' ]']); 
end

if success 
    theta = max(max(problem.proof.eps0,problem.proof.eps1),1);
    errorbound = altsup(max(theta.*rsol));
    disp('The proof was successful: there is a solution');
    if isfield(problem,'displaybounds') && problem.displaybounds
        disp(['The validation radii are ',num2str(rsol')]);
    end
    disp(['The uniform error bound is ',num2str(errorbound)]);    
    if altisintval(x(1))
        disp('Full proof including interval arithmetic')
    else
        disp('But this does not include interval arithmetic');
    end
else
    disp('The proof failed');
    rsol=NaN;
    errorbound=NaN;
end 

% storing some bounds
bounds.Y=altsup(Y);
bounds.Z=altsup(Z);
bounds.W=altsup(W);
bounds.eta=eta;

end


