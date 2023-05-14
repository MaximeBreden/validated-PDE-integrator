function y = altisintval(x)
% alternative to intlab-native isintval function
% to make it possible to run part of the code 
% without having access to intlab

global intervalarithmeticavailable 

if intervalarithmeticavailable 
    y=isintval(x);
else
    y=false;
end

end

