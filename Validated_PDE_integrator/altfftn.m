function y = altfftn(x)
% computes the FFT of a tensor of any dimensions
% for variables of double or intval type
%
% for intval type the size of the input (in all dimensions) 
% have to be powers of 2

if altisintval(x(1))
  n = size(x);
  m = length(n);
  y = x;   
  for dimension = 1:m
      if n(dimension)>1   % no fft if first dimension is singleton! 
          y = reshape(verifyfft(y,1),size(y));
      end
      y = permute(y,[2:m 1]);
  end
else
  y=fftn(x); 
end

return

