function Ginv = estimateGinv(x,problem)
% computes the bounds in the estimate of tail of the inverse of G

N=(size(x,1)-1)/2;
D=size(x,3);

Ginv=altzeros([D,D],x(1));
for d=1:D
    Ginv(d,d)=1;
    for l=1:d-1
        Ginv(d,l)=exp(boundmu(d,l,N,problem));
    end
end

end