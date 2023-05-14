function values=interpolationconstant(K,l,intvaltest)
% Computes the interpolation constant called \sigma_{K,l} in the paper, so
% that || f-P_K(f) || \leq \sigma_{K,l} || f^{(l+1)} ||,
% (with interval arithmetic if intvaltest is an interval).
% When l=0 and K<=6, gives the new constant obtained in the appendix of the
% paper (Theorem B.6). Otherwise, uses the formula given in Theorem B.1.
% 
% works for a vector of l values

llength=length(l);

if exist('intvaltest','var') && exist('intval','file') && isintval(intvaltest(1))
    useintervals=true;
    values=intval(zeros([llength,1]));
else
    useintervals=false;
    values=zeros([llength,1]);
end

if useintervals
    % precompute some factorials and related factors
    maxl=max(l);
    ifac=intval(zeros([maxl+1,1]));
    for i=0:maxl
        ifac(i+1)=ifactorial(i);
    end
    ifac2=ifac.^2;
    powfour=4.^intval((0:maxl)');
    pow4fac2=powfour.*ifac2;
    lebesgueK=lebesgueconstant(K,intval(1));
end

for lcount=1:llength
    ll=l(lcount);

    if ll > K
        disp('invalid inputs');
        val = NaN;
    elseif useintervals
        if ll==K 
            iK = intval(K);
            val = 1/(ifactorial(K+1)*2^(iK-1));
        elseif ll == 0
            switch K
                case 1
                    val = intval(1);
                case 2
                    val = (14*sqrt(intval(7))-20)/27;
                case 3
                    val = intval('0.6667');
                case 4
                    val = intval('0.5254');
                case 5
                    val = intval('0.4945');
                case 6
                    val = intval('0.4265');
                otherwise
                    ipi = intval('pi');
                    val = ipi/(2*(K+1))*(1+lebesgueK);
            end
        else     
            ipi = intval('pi');
            val1 = (ipi/2)^(ll+1)/prod(intval(K+1-ll:K+1))*(1+lebesgueK);
            val2 = 0;
            for q = 0:floor(ll/2)
                val2 = val2 + 1/(pow4fac2(q+1)*ifac(ll-2*q+1));
            end
            val2 = val2/(ll+1);
            val = min(val1,val2);
        end
    else
        if ll == K
           val = 1/(factorial(K+1)*2^(K-1));
        elseif ll == 0
            switch K
                case 1
                    val = 1;
                case 2
                    val = (14*sqrt(7)-20)/27;
                case 3
                    val = 0.6667;
                case 4
                    val = 0.5254;
                case 5
                    val = 0.4945;
                case 6
                    val = 0.4265;
                otherwise
                    val = pi/(2*(K+1))*(1+lebesgueconstant(K));
            end
        else
            val1 = (pi/2)^(ll+1)/prod(K+1-ll:K+1)*(1+lebesgueconstant(K));
            val2 = 0;
            for q = 0:floor(ll/2)
                val2 = val2 + 1/(4^q*factorial(ll-2*q)*factorial(q)^2);
            end
            val2 = val2/(ll+1);
            val = min(val1,val2);
        end
    end
    values(lcount)=val;
end

end % of function interpolationconstant
    
function Lambda = lebesgueconstant(K,intvaltest)
% Computes the Lebesque constant (with interval arithmetic if intvaltest is an interval) 
% When K is odd, we use an exact formula. When K is even, we use an
% upper-bound, except for small values (namely K = 2 and K = 4), for which
% we computed the exact value

if exist('intvaltest','var') && exist('intval','file') && isintval(intvaltest(1))
    ipi = intval('pi');
    i1 = intval(1);
    i2 = intval(2);
    sqrt2 = sqrt(i2);
    % computation needed for the case K = 4
    x = infsup(inf(intval('0.37227692859622')),sup(intval('0.37227692859624')));
%     x = verifynlss('8*x^3-6*(1+sqrt(0*x+2))*x^2-6*x+1+2*sqrt(0*x+2)',0.38);    
else
    sqrt2 = sqrt(2);
    ipi = pi;
    i1 = 1;
    i2 = 2;
    % computation needed for the case K = 4
    x = 0.37227692859623;
%     x = fzero('8*x^3-6*(1+sqrt(2))*x^2-6*x+1+2*sqrt(2)',0.38);
end

if K == 2
    Lambda = 5/i2^2;
elseif K == 4
    Lambda = 2*x^4-2*(1+sqrt2)*x^3-3*x^2+(1+2*sqrt2)*x+1;
else
    Lambda = sum(cot((1:2:2*K-1)*ipi/(4*K)))/K; %sharp if K odd, possible overestimation by at most 1/K^2 if K even
    % Lambda = 1+2/ipi*log(K+i1); % valid upperbound, but less sharp than the one above
end
    
end % of function lebesgueconstant
    
function val = ifactorial(n)
% interval arithmetic factorial 
if n == 0
    val = intval(1);
else
    val = prod(intval(1:n));
end

end % of function ifactorial 

