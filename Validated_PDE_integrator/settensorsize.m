function cout = settensorsize(cin,N)
% sets the size of the tensor to N
% for the output size(cout)==2*N+1
% It pads with zeros in dimensions that are too small
% truncates in dimensions that are too large 

dim=length(N);
M = (size(cin)-1)/2;
M(length(M)+1:dim) = 0;

L = max(N,M)+1;

clarge = altzeros(2*L-1,cin(1));
sin=cell(dim,1);
sout=cell(dim,1);
for j=1:dim
    sin{j}=L(j)+(-M(j):M(j));
    sout{j}=L(j)+(-N(j):N(j));
end

clarge(sin{1:dim}) = cin;
cout = clarge(sout{1:dim});          

if dim==1
% turn row into column
    cout=cout.';
end

end


