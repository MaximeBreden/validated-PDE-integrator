function out = obj_func_fl(X,K,M)

% The objective function to be minimized, for the first (nonrigorous) step
% where we look for a candidate minimizer. P corresponds to P_K(y)(t) in 
% the paper, and w to the value of the function being interpoated at t.

t = X(1,:); 
w = X(2,:);
y = X(3:end,:);

P = M*y;
TK = cos((0:K)' * acos(t));
P = sum(P.*TK,1);
out = - abs(P-w);