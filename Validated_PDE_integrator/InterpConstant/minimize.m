function [current_min, minimizers] = minimize(K,M,x,dx,box,nb_sub,tol,depth,depth_max,threshold,lsub,lmax,keeptrackofminimizers,show)

% Function which will be called recursively in order to enclose the
% minimum, and to find potential minimizers (more precisely, boxes whose 
% reunion must containt all minimizers).

if show
    fprintf("\nCurrent depth: %d\n",depth)
end
minimizers = [];

%% Subdivision of the input into smaller boxes
X = subdivide_adapt(box,K,nb_sub,lmax); 
if show
    fprintf("Numbers of elements after subdivison : %d\n",size(X,2))
end

%% Keeping only the columns of X compatible with ||f'||<=1
X = check_input(X,dx); 
if show
    fprintf("Numbers of admissible elements : %d\n",size(X,2))
end
if isempty(X)
    current_min = Inf;
    return
end

%% We evaluate the objective function on each column of X
fX = obj_func(X,K,M,x);

%% Sort the columns of X according to the values of f(X), and discard those which cannot contain a minimizer
[X, fX, current_min] = check_output(X,fX,threshold); 
if show
    fprintf("Numbers of elements in which the minimum could be : %d\n",size(X,2))
    infsup(current_min)
end

if isempty(X)
    return
end

%% Calling the minimize function again on X
threshold = min(threshold,sup(current_min));
if rad(current_min)>tol && depth<depth_max
    L = size(X,2);
    if L <= lsub % X is not too large, we reapply minimze to the whole X
        [current_min, minimizers] = minimize(K,M,x,dx,X,nb_sub,tol,depth+1,depth_max,threshold,lsub,lmax,keeptrackofminimizers,show);
    else % X is too large (with respect to lsub), we split X into blocks of lsub rows and reapply minimize to each of them
        J = ceil(L/lsub);
        tab_min_loc = intval(Inf*ones(1,J));
        j = 1;
        Xj = X(:,(j-1)*lsub+(1:lsub)); %first block
        fXj = min(fX((j-1)*lsub+(1:lsub)));
        threshold_start = threshold;
        while j<=ceil(L/lsub) && not(threshold<fXj)
            [current_min_loc, minimizers_loc] = minimize(K,M,x,dx,Xj,nb_sub,tol,depth+1,depth_max,threshold,lsub,lmax,keeptrackofminimizers,show);
            tab_min_loc(j) = current_min_loc;
            minimizers = [minimizers minimizers_loc];
            threshold = min(threshold,sup(current_min_loc));
            j = j+1;
            if j<J
                Xj = X(:,(j-1)*lsub+(1:lsub)); %next block
                fXj = min(fX((j-1)*lsub+(1:lsub)));
            else %the last block might be smaller
                Xj = X(:,(j-1)*lsub+1:end); 
                fXj = min(fX((j-1)*lsub+1:end));
            end         
        end  
        current_min = min(tab_min_loc);
        if keeptrackofminimizers && threshold < threshold_start 
            minimizers = check_output(minimizers,obj_func(minimizers,K,M,x),threshold); 
        end
    end
elseif keeptrackofminimizers
    minimizers = X;
end
