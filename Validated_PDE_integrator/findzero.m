function xout = findzero(xin,problem)
% performs Newton iterations to find a zero of the problem F=0 starting at xin

x=xin;
N=(size(x,1)-1)/2;
K=size(x,2)-1;
D=size(x,3);
dim=(2*N+1)*(K+1)*D;

problem.type='numerics';

% initialize Newton
nsteps=24; 
tolerance=1e-13;
normdx=ones(nsteps+1,1);
n=1;

% delay is used to check it does not aimlessly wander
% with small numbers larger than tolerance but smaller than relaxedtol
delay=2;
relaxedtol=1e-11;
% for default behaviour choose delay=0; relaxedtol=tolerance;
% delay=0;
% relaxedtol=tolerance;

disp('Start of Newton iterations')
while n<=nsteps && normdx(n)>=tolerance && ~(n>delay && all(normdx(n-delay:n)<relaxedtol))
    F=determineF(x,problem);
    DF=determineDF(x,problem);
    F=reshape(F,dim,1);
    DF=reshape(DF,dim,dim);
    dx=-DF\F;
    % enforce symmetry
    if isfield(problem,'symmetry') 
        if strcmp(problem.symmetry,'cosineseries')
            dx=real(dx);
        elseif strcmp(problem.symmetry,'sineseries')
            dx=1i*imag(dx);
        end
    end
    n=n+1;
    normdx(n)=norm(dx);
    disp(['stepsize is ',num2str(normdx(n))]);
    x=x+reshape(dx,[2*N+1,K+1,D]);
end

xout=x;
disp('End of Newton iterations')

end

