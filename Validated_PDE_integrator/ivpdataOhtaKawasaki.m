function [problem,plotdata] = ivpdataOhtaKawasaki
% data for Ohta-Kawasaki equation
%
% u_t = - 1/gamma^2 u_{xxxx} - (u-u^3)_{xx} - sigma(u-m)
%
% for x in [0,L] and t in [0,tau]
%
% returns the nonlinearities, 
% the problem definition, including the initial data 
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
% u_t = - 1/gamma^2 u_{xxxx} - (u-u^3)_{xx} - sigma(u-m)
% for x in [0,L] and t in [0,tau]

gamma=sqrt(8*altone); 
sigma=altone/5;
m=altone/10;
L=2*altpi;
tau=5*altone;

% rescale to spatial domain [0,2*pi] and
% u_t = -u_{xxxx} + lambda21*(u-u^3)_{xx} + lambda01*(u-m)

scale.space=L/(2*altpi);
scale.time=gamma^2*scale.space^4;
integrationtime=tau/scale.time;

lambda21=-scale.time/scale.space^2;
lambda23=-lambda21;
lambda01=-sigma*scale.time;
lambda00=-lambda01*m;

%% Truncation and grid %%

% number of Fourier modes (2N+1 really)
N=21; 
% number of Chebyshev modes (K+1 really)
K=2;
% number of time domains (called M in the latex)
D=395; 

% total integration time
griddata.T=integrationtime;
griddata.D=D;
% grid adjustment with linear (skew) and quadratic (dip) term
% skew=dip=0 means no skew (uniform grid)
% skew>0 means smaller steps at start
% dip>0 means smaller steps in the middle
% range is roughly -1<skew<1 and -2<dip<4
griddata.skew=-0.7;  
griddata.dip=0.8; 

%% Initial data %%

% set initial data 
initialdata=altzeros([2*N+1,1],m);
initialdata(N+1)=m;
initialdata(N+2)=altone/10;

% note that these will be symmetrized later as follows (note factor 2):
% initialdata=(initialdata+conj(initialdata(end:-1:1)))/2;

%% Parameters for the proof %%

% weigths in the norms for the Banach space
nu=1.0001;
eps0=3; % eps_{infty N} in the latex
eps1=0.4; % eps_{infty} in the latex
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
  alllambda=[lambda21,lambda23,lambda01,lambda00];
  epsfactor=max(sup(intval(rad(alllambda))./eps(mid(alllambda))));
  % turn relevants coefficient(s) into floats
  lambda21=mid(lambda21);
  lambda01=mid(lambda01);
  m=mid(m);
end

% nonlinearity{j} represented the term g^{(j-1)}(u)
nonlinearity{1} = lambda01*(u-m);
nonlinearity{3} = lambda21*(u-u^3);

%% Preprocessing %%
[problem.pde.polynomials,symmetry]=preprocessnonlinearities(nonlinearity,epsfactor);
problem = preprocessdata(problem,initialdata,griddata,symmetry,scale);
problem.K = K;

%% Final settings %%
problem.useintervals = [];

%% Plot data %%
plotdata=[];


end