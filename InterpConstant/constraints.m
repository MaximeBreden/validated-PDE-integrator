function [c,ceq] = constraints(X,K,x,dx,tol)

% The constraints, for the first (nonrigorous) step where we look for a 
% candidate minimizer. w corresponds to the value of the function f being 
% interpolated at t, and y = (y_0,...,y_K) to the values taken by f at the 
% interpolation nodes x = (x_0,...,x_K). The constraints all stem from the 
% fact that f(0) = 0 and ||f'|| <= 1, which implies that 
% |y_{k+1}-y_k| <= |x_{k+1}-x_k|, |w| <= |t|, and |w-y_k| <= |t-x_k|. This
% last constraint only needs to be enforced for half the k's, since by
% symmetry we restrict t from [-1,1] to [-1,0].

t = X(1,:);
w = X(2,:);
y = X(3:end,:);

dy = abs(y(2:end,:)-y(1:end-1,:));
ind = ceil(K/2);

c = [dy-dx; % |y_{k+1}-y_k| - |x_{k+1}-x_k|
     abs(w)-abs(t); % |w| - |t|
     abs(w-y(1:ind,:))-abs(t-x(1:ind))]; % |w-y_k| - |t-x_k| (for an appropriate subset of indices)
if nargin == 5
    c = c - tol;
end
ceq = [];