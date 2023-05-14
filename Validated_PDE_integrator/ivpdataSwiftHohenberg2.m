function [problem,plotdata] = ivpdataSwiftHohenberg2
% data for Swift-Hohenberg equation
%
% u_t = - u_{xxxx} - 2 u_{xx} + (alpha-1)*u - u^3
%
% for x in [0,L] and t in [0,tau]
%
% returns the problem definition, including the initial data 
% with the initial data setting the number N of Fourier modes
% the number K of Chebyshev modes (interpolation nodes) in time
% the number D of domain in the (time) domain decomposition
% (note that D is called M in the latex)
%
% the problem definition also includes some constants needed in the proof
%
% additional data for the plot can also be specified

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
% u_t = - u_{xxxx} - 2 u_{xx} - u + alpha*u - u^3
% for x in [0,L] and t in [0,tau]

alpha=5*altone; 
L=6*altpi;
tau=3*altone/2;

% rescale to spatial domain [0,2*pi] and
% u_t = -u_{xxxx} +lambda2*u_{xx} + lambda01*u+lambda03*u^3

scale.space=L/(2*altpi);
scale.time=scale.space^4;
integrationtime=tau/scale.time;

lambda2=-2*scale.time/scale.space^2;
lambda01=(alpha-1)*scale.time;
lambda03=-scale.time;

%% Truncation and grid %%

% number of Fourier modes (2N+1 really)
N=19; %N=20
% number of Chebyshev modes (K+1 really)
K=2;
% number of time domains (called M in the latex)
D=108; %D=103

% total integration time
griddata.T=integrationtime;
griddata.D=D;
% skew=0 means no skew (uniform grid)
% skew>0 means smaller steps at start, range is between -1 and 1
griddata.skew=0.4;  

%% Initial data %%

% set initial data 
initialdata=zeros(2*N+1,1);
initialdata(N+2)=0.4;
initialdata(N+3)=-0.3;

% note that these will be symmetrized later as follows (note factor 2):
% initialdata=(initialdata+conj(initialdata(end:-1:1)))/2;

%% Parameters for the proof %%

% weigths in the norms for the Banach space
nu=1.0001;
eps0=0.2; % eps_{infty N} in the latex
eps1=0.1; % eps_{infty} in the latex
% nu is the same on all domains 
problem.proof.nu=nu;
% in principle the epsilon weights may differ per domain
problem.proof.eps0=repmat(eps0,D,1); 
problem.proof.eps1=repmat(eps1,D,1);

% a apriori bound on the validation radius
rstar=2e-1;
% this may also be different on each domain
problem.proof.rstar=repmat(rstar,D,1);

% computational constants used in computing exponential integrals:
% the number of quadrature nodes for computing numerically 
problem.numerics.integrals.K1=100; 
% the number of quadrature nodes used for comutations in the proof
problem.proof.integrals.K1=100; 
% the number of Chebyshev modes used to represent the integrand in the proof
problem.proof.interpolation.K0=50; 

%% PDE definition %%

% symbolic variable
syms u;

% the (even) order of the PDE (2R in the latex)
problem.pde.order=4;

if intervalarithmeticavailable
  % determine maximum error in rounding, relative to eps
  % for all monomial coefficients
  alllambda=[lambda2,lambda01,lambda03];
  epsfactor=max(sup(intval(rad(alllambda))./eps(mid(alllambda))));
  % turn coefficient(s) into floats
  lambda2=mid(lambda2);
  lambda01=mid(lambda01);
  lambda03=mid(lambda03);
end

% nonlinearity{j} represented the term g^{(j-1)}(u)
nonlinearity{1} = lambda01*u+lambda03*u^3;
nonlinearity{3} = lambda2*u;

%% Preprocessing %%
[problem.pde.polynomials,symmetry]=preprocessnonlinearities(nonlinearity,epsfactor);
problem = preprocessdata(problem,initialdata,griddata,symmetry,scale);
problem.K = K;

%% Final settings %%
problem.useintervals=[];

%% Plot data %%
plotdata=[];

end