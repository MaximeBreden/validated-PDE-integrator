function matrix = cheb_evalmatrix(t,K)
% matrix corresponding to evaluating the Chebyshev polynomial 

k=(0:K)';
matrix=2*cos(k*acos(t));
matrix(1,:)=1;

end

