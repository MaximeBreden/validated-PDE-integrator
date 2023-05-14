function x=numericalintegration(problem)
% Numerical integration of the Fourier-truncated PDE
% to find a guess for the solution
% The number of Chebyshev nodes is K+1
% The initial data are stored in problem.initial.data

initial=problem.initial.data;
N=(length(initial)-1)/2;
D=length(problem.timegrid)-1;
K=problem.K;
domains=diff(problem.timegrid);
x=zeros(2*N+1,K+1,D);

% no intervals; this is numerics
for p=1:length(problem.pde.polynomials)
    problem.pde.polynomials{p}=altmid(problem.pde.polynomials{p});
end

for d=1:D
    % in each domain we want the solution at the Chebyshev nodes
    % (which are ordered from 1 to -1, hence the two "fliplr")
    times=(fliplr(cheb_grid(K))+1)/2*domains(d);
    [~,orbit]=ode45(@(t,a) vectorfield(a,problem),times,initial);
    if K==1
        % extracting begin and end point
        orbit=orbit([1 end],:);
    end
    x(:,:,d)=cheb_coeffs(fliplr(orbit.'));
    % initial data for the next domain
    initial=orbit(end,:).';
end

end

function f=vectorfield(a,problem)
% vectorfield of the ODE which results from truncating the PDE
    N=(length(a)-1)/2;
    % leading term
    f=(-((-N:N)').^problem.pde.order).*a;
    for p=0:length(problem.pde.polynomials)-1
        % terms with p-th derivative
        cp=convnonlinearity(a,p,0,problem);
        cp=diag((1i.*(-N:N)).^p)*cp;
        f=f+cp;
    end
end
