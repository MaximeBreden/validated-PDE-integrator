function [aout1,aout2] = convnonlinearity(ain,k,m,problem)
% evaluates convolution nonlinearity of
% D_u^m g_k(u) for m=0,1 (which get k derivatives of x in the PDE)
% and the termwise absolute value of D_u^2 g_k(u) for m=2

[aout1,aout2] = convpolynomialtensor(ain,@(u) nonlinpoly(u,k,m,problem),nonlindegree(k,m,problem));

end
