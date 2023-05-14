function z = nonlinpoly(y,k,m,problem,absolute)
% evaluate nonlinearity corresponding to
% D_u^m g_k(u) for m=0,1 (which get k derivatives of x in the PDE)
% and the termwise absolute value of D_u^2 g_k(u) for m=2
% k will never exceed (or be equal to) problem.pde.order

if k+1>length(problem.pde.polynomials) || isempty(problem.pde.polynomials{k+1})
    % in case the term is not defined it vanishes
    polynomial=[];
else
    if altisintval(y(1))
        polynomial=problem.pde.polynomials{k+1};
    else
        % no need to compute with interval coefficients
        polynomial=altmid(problem.pde.polynomials{k+1});
    end
end

if exist('absolute','var') && strcmp(absolute,'absolute')
    % absolute value of coefficients
    polynomial=abs(polynomial);
end

if m==0
    z=evalpoly(polynomial,y);
elseif m==1
    n=length(polynomial)-1;
    if n>0
        Dp=polynomial.*(n:-1:0);
    else
        % catching some problem with older matlab versions
        Dp=[];
    end        
    z=evalpoly(Dp(1:n),y);
elseif m==2
    n=length(polynomial)-1;
    if n>1
        DDp=polynomial.*(n:-1:0).*(n-1:-1:-1);
    else
        DDp=[];
    end        
    z=evalpoly(DDp(1:n-1),y);
else
    disp('m should not be larger than 2');
    z=y*NaN;
end

end

function y=evalpoly(p,x)
% evaluates polynomial p in points x
  if isempty(p)
      p=altzeros([1,1],x);
  end
  y=repmat(p(end),size(x));
  for j=1:length(p)-1
      y=y+p(end-j).*x.^j;
  end
end
