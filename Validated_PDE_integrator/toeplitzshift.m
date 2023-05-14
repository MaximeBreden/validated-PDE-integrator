function matrix = toeplitzshift(v)
% toeplitz matrix with v in the middle column

N=(length(v)-1)/2;
column=[v(N+1:2*N+1);altzeros([N,1],v(1))];
row=[v(N+1:-1:1).',altzeros([1,N],v(1))];
matrix=toeplitz(column,row);

end

