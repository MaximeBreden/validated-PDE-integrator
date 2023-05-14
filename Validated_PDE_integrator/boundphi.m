function [phi,checkphi,tildephi] = boundphi(a,b,N,eps0,eps1,nu)
% computes the convolution bounds on tail elements 
%
% phi and checkphi are columns of length 2*N+1, tildephi is a scalar
%
% the input vectors a and b are bounds on the elementswise C_0 norms

Na = (length(a)-1)/2;
Nb = (length(b)-1)/2;

phi = altzeros([2*N+1,1],a(1));
checkphi = altzeros([2*N+1,1],a(1));

invnu=(1/nu).^(0:N+Nb)';

for n = -N:N
   ma = max(-Na,n-N):min(Na,n+N);
   amax = max(a(Na+1+ma).*invnu(1+abs(n-ma)));

   mb = [-Nb:n-N-1,n+N+1:Nb]; 
   if isempty(mb)
        bmax=0;
   else
        bmax = max(b(Na+1+mb).*invnu(1+abs(n-mb)));
   end
   
   % take the max
   phi(N+1+n) = max(eps0*amax,eps1*bmax); 
   checkphi(N+1+n) = max(max(1,eps0)*amax,eps1*bmax); 
end

asummax = 0;
bsummax = 0;
for m = -N:N
   % n should be in the range where shifted a is nonzero/defined -Na+m:Na+m
   % and where we need to take the sum, namely -N:N    
   n = max(-N,-Na+m):min(N,Na+m);
   asum = sum(a(Na+1+n-m).*nu.^(abs(n)'-abs(m)));
   asummax = max(asummax,asum);     
end
% m with |m|>N should be in the range where we will encounter a nonzero b
% namely between -N-Nb and N+Nb
for m = [-N-Nb:-N-1,N+1:N+Nb]
   % n should be in the range where shifted a is nonzero/defined -Nb+m:Nb+m
   % and where we need to take the sum, namely -N:N    
   n = max(-N,-Nb+m):min(N,Nb+m);
   bsum =  sum(b(Nb+1+n-m).*nu.^(abs(n)'-abs(m)));
   bsummax = max(bsummax,bsum);     
end
% take the final max
tildephi = max(max(1,eps0)*asummax,eps1*bsummax);

end

