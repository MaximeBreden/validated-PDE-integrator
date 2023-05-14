function y=altmid(x)
% takes mid of interval input
% does nothing for double input

if isempty(x) 
    y=[];
elseif altisintval(x(1))
    y=mid(x);
else
    y=x;
end
