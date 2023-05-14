function [problem,plotdata] = ivpdataFisher3
% data for Fisher's equation
% with naive semigroup
%
% u_t = mu*u_{xx} + lambda*u*(1-u)
%
% for x in [0,L] and t in [0,tau]
%
% returns the problem definition, including the initial data 
% with the initial data setting the number N of Fourier modes
% the number K of Chebyshev modes (interpolation nodes) in time
% the number D of domains in the (time) domain decomposition
% (note that D is called M in the latex)
%
% the problem definition also includes some constants needed in the proof
%
% additional data for the plot can also be specified

%% !! Naive Semigroup !!
problem.semigroup.type='naive';

global intervalarithmeticavailable 

if intervalarithmeticavailable
    altpi=intval('pi');
    altone=intval(1);
else
    altpi=pi;
    altone=1;
    epsfactor=[];
end

%% Set problem parameters %%

% define the constants in the equation
% u_t = mu*u_{xx} + lambda*u*(1-u)
% for x in [0,L] and t in [0,tau]

mu=altone; % diffusion constant
lambda=altone;
L=4*altpi;
tau=4*altone;
    
% rescale to spatial domain [0,2*pi] and
% u_t = u_{xx} + lambda0*u*(1-u)

scale.space=L/(2*altpi);
scale.time=scale.space^2/mu;
integrationtime=tau/scale.time;

lambda0=lambda*scale.time;

%% Truncation and grid %%

% number of Fourier modes (2N+1 really)
N=14;
% number of Chebyshev modes (K+1 really)
K=2;
% number of time domains (called M in the latex)
D=165; 

% total integration time
griddata.T=integrationtime;
griddata.D=D;
% =0 means no skew (uniform grid); =1 is maximum skew
griddata.skew=0.9;  

%% Initial data %%

% set initial data 
initialdata=zeros(2*N+1,1);
initialdata(N+1)=0.5;
initialdata(N+2)=1i*0.5;
initialdata(N+3)=-1;
initialdata(N+5)=-1i*0.2;

% note that these will be symmetrized later as follows (note factor 2):
% initialdata=(initialdata+conj(initialdata(end:-1:1)))/2;

%% Parameters for the proof %%

% weigths in the norms for the Banach space
nu=1.0001;
eps0=0.3; % eps_{infty N} in the latex
eps1=0.3; % eps_{infty} in the latex
% nu is the same on all domains 
problem.proof.nu=nu;
% in principle the epsilon weights may differ per domain
problem.proof.eps0=repmat(eps0,D,1); 
problem.proof.eps1=repmat(eps1,D,1);

% a apriori bound on the validation radius
rstar=7e-1;
% this may also be different on each domain
problem.proof.rstar=repmat(rstar,D,1);

% computational constants used in computing exponential integrals:
% the number of quadrature nodes for computing numerically 
problem.numerics.integrals.K1=50; 
% the number of quadrature nodes used for comutations in the proof
problem.proof.integrals.K1=50; 
% the number of Chebyshev modes used to represent the integrand in the proof
problem.proof.interpolation.K0=50; 
% note that in the latex K0=K_1 and K1=K_0

%% PDE definition %%

% symbolic variable
syms u;

% the (even) order of the PDE (2R in the latex)
problem.pde.order=2;

if intervalarithmeticavailable
  % determine maximum error in rounding, relative to eps
  epsfactor=sup(intval(rad(lambda0))/eps(mid(lambda0)));
  % turn coefficient(s) into floats
  lambda0=mid(lambda0);
end

% nonlinearity{j} represented the term g^{(j-1)}(u)
nonlinearity{1} = lambda0*(u-u^2);

%% Preprocessing %%
[problem.pde.polynomials,symmetry]=preprocessnonlinearities(nonlinearity,epsfactor);
problem = preprocessdata(problem,initialdata,griddata,symmetry,scale);
problem.K = K;

%% Final settings %%
problem.useintervals = [];

%% Plot data %%
plotdata=[];

end
