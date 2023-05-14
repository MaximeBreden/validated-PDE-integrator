clear variables
close all
clc

% This script produces the proof of Theorem B.6. You just have to change
% the value of K below in order to get the different cases.

% Convention: since the code is written for minimization problems, we
% replace max_{[-1,1]\times\Omega_K} max( \overline{\phi}(y)(t)-P(y)(t), P(y)(t) - \underline{\phi}(y)(t) )
% by - min_{[-1,1]\times\Omega_K} min( P(y)(t)-\overline{\phi}(y)(t), \underline{\phi}(y)(t) - P(y)(t) )
% and solve the latter.

% /!\ You cannot run this script without Intlab /!\

%% Initialization
K = 2; %The K for which we want to compute \sigma_{K,0}
nb_opt_nonrig = 0; %Number of random initial data that are tried by the nonrigorous solver

[box,threshold,M,x,dx] = initialize_optim(K,nb_opt_nonrig);

%% Optimization by "bissection"

tic

% Parameters
keeptrackofminimizers = false; % Put to "true" if you want to keep track of all the boxes in which a minimizer could be (can create memory issues and significantly increase runtime for larger values of K)
nb_sub = 6; % Number of pieces in which the boxes are split (component-wise) when subdividing
if K <= 6
    depth_max = 14; % The input will be subdivided recursively up to depth_max times
else
    depth_max = 12; % The input will be subdivided recursively up to depth_max times
end
% depth_max is one of the critical parameters influencing the precision of
% the output, but also the time taken by the algorithm. To be adjusted
% depending on your needs an the computational power at your disposal.
tol = 1e-5; % If the interval enclosing the minimal value has a radius <= tol, we stop
lmax = 2e4; % Controls the numbers of elements that can be produced by a subdivision step (overriding nb_sub if needed)
lsub = 2e4; % If there are more that lsub elements in a subdivision, we do not deal with all of them at once

show = 1; % Controls the intermediate outputs 

% If you already know an upperbound for the minimal value (for instance
% from a previous computation with smaller depth_max), or from some pen and
% paper computations, you can input this threshold below in 
% "known_threshold", which is put to Inf by default.
known_threshold = Inf;
threshold = min(threshold,known_threshold);

fprintf("\nStarting the rigorous minimization\n")
depth_ini = 0; 
% This is where things happen
[mu,minimizers] = minimize(K,M,x,dx,box,nb_sub,tol,depth_ini,depth_max,threshold,lsub,lmax,keeptrackofminimizers,show); 
mu = infsup(inf(mu), min(sup(mu),threshold));
mu = -mu;
fprintf("\nThe constant sigma_{K,0} for K = %d is contained in\n",K)
infsup(mu)

toc

%% Potentially improving the data stored for a future run
if not(isempty(minimizers))
    X = minimizers(:,1);
    load('data_opt.mat','Xstored')
    str = ['K',num2str(K)];
    if isfield(Xstored,str)
        XX = Xstored.(str);
        if sup(obj_func(X,K,M,x)) < sup(obj_func(XX,K,M,x))
            Xstored.(str) = X;
            save('data_opt.mat','Xstored')
        end
    else
        Xstored.(str) = X;
        save('data_opt.mat','Xstored')
    end
end



