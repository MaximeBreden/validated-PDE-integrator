function z = altzeros(n,intvaltest)
% Z = altzeros(N,INTVALTEST)
% produces a tensor of zeros of double or interval type
% depending on the type of INTVALTEST
% You must (really: must) bracket the dimensions in N, i.e. 
% altzeros([3,4],intval('3.1'))
% The default is no intervals (Z is of type double if INTVALTEST is not provided)

if length(n)==1
  disp('please check code; possible problem in altzeros; you have to bracket dimensions')
end

z = zeros(n);

if exist('intvaltest','var') && (altisintval(intvaltest(1)) || ...
                                (islogical(intvaltest) && intvaltest))
  z=intval(z);
end

return

