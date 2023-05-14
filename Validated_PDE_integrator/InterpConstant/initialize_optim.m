function [box,threshold,M,x,dx] = initialize_optim(K,nb_opt_nonrig)

% Computes box = { y in R^{K+1}, |y_k| <= |t_k| for k=0..K} on which we 
% have to minimize the objective function. The additional constraints 
% |y_{k+1}-y_{k}| <= |t_{k+1}-t_{k}| are taken care of within the function
% check_input.m
% Also runs a non-rigorous optimization (on a slightly different version of
% the problem), with nb_opt_nonrig randomly generated initials points, in
% order to find threshold, which is an upperbound on the minimal value of 
% the objective function.

ind = 0:K;
theta = ((K-ind')*intval(pi))/K;
x = cos(theta); % the Chebyshev nodes

% Construction of the matrix to go from values at Chebyshev nodes to
% Chebyshev coefficients
M = transpose( cos( theta * ind ) ) / K; 
M(:,[1,K+1]) = M(:,[1,K+1]) / 2;
M(2:K,:) = 2 * M(2:K,:);

bounds_coeffs = sup(abs(x));
dx = abs(x(2:end,:)-x(1:end-1,:));

%% Creating the box on which to optimize
box = intval(zeros(K+2,1));
box(1) = infsup(-1,0); % bounds for t
box(2) = infsup(0,bounds_coeffs(1)); % we can always assume f(-1) >= 0
for k = 2 : K+1
    box(k+1) = midrad(0,bounds_coeffs(k)); % bounds for the value f(x_k) (namely |f(x_k)| <= |x_k|)
end

%% Getting a first upper bound on the minimal value, either by loading data from previous runs or by using nonrigorous optimization

threshold = Inf;
str = ['K',num2str(K)];

if exist('data_opt.mat','file') == 2
    load('data_opt.mat','Xstored')
    if isfield(Xstored,str)
        fprintf('\nStored potential minimizer: ')
        X = Xstored.(str)
        fprintf('\nThe minimal value cannot be larger than: ')
        threshold = sup(obj_func(X,K,M,x))
    end
end

if exist('nb_opt_nonrig','var') && nb_opt_nonrig > 0
    fprintf("\nNon rigorous resolution of the optimization problem, used to get a rigorous upper bound on the minimal value\n")
    % For finding a candidate for the minimizers, we use a slightly
    % reformulated version of the problem, making it easier to use Matlab's
    % constrained nonlinear optimizer fmincon. That is, rather than
    % computing the upper and lower bounds (\overline\varphi and
    % \underline\varphi in the paper), we add an extra variable w
    % representing the value of the function to be interpolated, at the
    % point t, and also optimize over w, with added constraints ensuring
    % that the value w is compatible with the other values taken at the
    % interpolation nodes, and the fact that ||f'|| <= 1.
    Xopt = [];
    valopt = Inf;
    ddx = inf(dx);
    box2 = [box(1);infsup(-1,1);box(2:end)];    
    f = @(X) obj_func_fl(X,K,mid(M));
    c = @(X) constraints(X,K,mid(x),mid(dx),2*eps);
    options1 = optimoptions(@fmincon,'Display','off');
    options2 = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
    options3 = optimoptions(@fmincon,'Algorithm','active-set','Display','off');
    warning('off','optimlib:fwdFinDiffInsideBnds:StepReduced')
    for j = 1:nb_opt_nonrig
        fprintf("Try %d out of %d ... ",j,nb_opt_nonrig)
        valid_Xin = false;
        it = 0;
        Xin = zeros(K+2,1);
        while not(valid_Xin) && it<1e6
            Xin(1) = -rand(1); %t in [-1,0]
            Xin(2) = abs(Xin(1)) * (2*rand(1)-1); %|w| <= |t|
            Xin(3) = rand(1); %y_0 in [0,1]
            for k = 2:K+1
                low = max(inf(box2(k+2)),Xin(k)-ddx(k-1));
                up = min(sup(box2(k+2)),Xin(k)+ddx(k-1));
                Xin(k+2) = low + (up-low)*rand(1);
            end
            valid_Xin = min(constraints(Xin,K,mid(x),mid(dx))<=0);
        end
        if valid_Xin
            fprintf("Valid initial data found\n")
            [Xopt1, val1] = fmincon(f,Xin,[],[],[],[],inf(box2),sup(box2),c,options1);
            [Xopt2, val2] = fmincon(f,Xin,[],[],[],[],inf(box2),sup(box2),c,options2);
            [Xopt3, val3] = fmincon(f,Xin,[],[],[],[],inf(box2),sup(box2),c,options3);
            if val1 < valopt
                valopt = val1;
                Xopt = Xopt1;
            end
            if val2 < valopt
                valopt = val2;
                Xopt = Xopt2;
            end
            if val3 < valopt
                valopt = val3;
                Xopt = Xopt3;
            end
        else
            fprintf("No valid initial data found\n")
        end
    end

    Xopt = Xopt([1,3:K+3]);
    rad = 1e-16;
    Xopt = midrad(Xopt,rad); %the point used to get the upper bound on the minimal value
    X = check_input(Xopt,dx);
    it = 1;
    % If needed, we enlarge Xopt to ensure that it contain some elements
    % satisfying the additional constraints |y_{k+1}-y_{k}| <= |t_{k+1}-t_{k}|
    while isempty(X) && it < 16
        rad = 10*rad;
        Xopt = midrad(mid(Xopt),rad);
        X = check_input(Xopt,dx);
        it = it+1;
    end

    if not(isempty(X))
        if sup(obj_func(X,K,M,x)) < threshold
            fprintf('\nBest minimizer found (nonrigorously): ')
            X
            fprintf('\nThe minimal value cannot be larger than: ')
            threshold = sup(obj_func(X,K,M,x))
            Xstored.(str) = X;
            save('data_opt.mat','Xstored')
        else
            fprintf('\nDid not find a better candidate than what was already stored\n')
        end
    end
end

