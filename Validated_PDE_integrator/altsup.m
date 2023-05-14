function y=altsup(x)
% takes sup of interval input
% does nothing for double input

if isempty(x) 
    y=[];
elseif altisintval(x(1))
    y=sup(x);
else
    y=x;
end
