function [y,problem] = ivpsolve(problem)
% solve the IVP numerically using an integrator 
% and then improve using Newton iterations
%
% The default is to use the more elaborate semigroup described in the paper.
% If problem.semigroup.type='naive' then 
% the problem setup uses the naive semigroup
% based on the highest order derivative only
% and ignoring any nonlinear terms.
%
% outputs symmetric data if the problem is symmetric

% numerical integation of IVP
disp('Numerical integration');
problemfloat=problem;
problemfloat.initial.data=altmid(problemfloat.initial.data);
problemfloat.timegrid=altmid(problemfloat.timegrid);
x=numericalintegration(problemfloat);

if isfield(problem,'symmetry') 
    if strcmp(problem.symmetry,'cosineseries')
        x=real(x);
    elseif strcmp(problem.symmetry,'sineseries')
        x=1i*imag(x);
    end
end

% Newton iteration with naive semigroup
disp('Solve with naive semigroup');
problemfloat.semigroup.type='naive';
problemfloat.semigroup = fixsemigroup(x,problemfloat);
y=findzero(x,problemfloat);

% Now repeat with elaborate semigroup (unless naive one is wanted)
if ~isfield(problem,'semigroup') || ~isequal(problem.semigroup.type,'naive') 
    residue=determineF(y,problemfloat);
    disp(['norm of the residue is ',num2str(norm(residue(:)))]);
    % Newton iteration with improved semigroup
    disp('Change the semigroup');
    problemfloat.semigroup.type='elaborate';
    problemfloat.semigroup = fixsemigroup(y,problemfloat);
    y=findzero(y,problemfloat);
end

% symmetrize
y=(y+conj(y(end:-1:1,:,:)))/2;
if isfield(problem,'symmetry') 
    if strcmp(problem.symmetry,'cosineseries')
        y=real(y);
    elseif strcmp(problem.symmetry,'sineseries')
        y=1i*imag(y);
    end
end

residue=determineF(y,problemfloat);
disp(['norm of the residue is ',num2str(norm(residue(:)))]);
problem.semigroup=problemfloat.semigroup;

end