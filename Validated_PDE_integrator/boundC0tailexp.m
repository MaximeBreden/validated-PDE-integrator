function bound = boundC0tailexp(N,d,problem)
% bound on exp(tau(d)*real(tlambda_n)*(t+1)) for all |n|>N and -1<=t<=1

%t=-1
bound0=1; 

%t=1
bound1=boundat1tailexp(N,d,problem);

%max attained at either t=-1 or t=1
bound=max(bound0,bound1);

end

