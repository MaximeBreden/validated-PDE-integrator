function y = expm1div(x)
% computes (exp(x)-1)/x
% when x is a float or an interval 
% improves on the standard algorithm for small x

if altisintval(x(1))
    x1=intval(inf(x));
    x2=intval(sup(x));
    w1=1+x1/2+x1.^2/6-abs(x1).^3/24./(1-abs(x1));
    w1(abs(x1)>=1)=-Inf;
    w2=1+x2/2+x2.^2/6+abs(x2).^3/24./(1-abs(x2));
    w2(abs(x2)>=1)=Inf;
    z1=(exp(x1)-1)./x1;
    z1(x1==0)=1;
    z2=(exp(x2)-1)./x2;
    z2(x2==0)=1;
    y1=max(inf(w1),inf(z1));
    y2=min(sup(w2),sup(z2));
    y=infsup(y1,y2);
else
    y=expm1(x)./x;
    y(x==0)=1;
end

