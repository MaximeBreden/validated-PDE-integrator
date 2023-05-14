function [fv,fvfull] = convpolynomialtensor(v,poly,order,parameters,zerocomponent)
% [FV,FVFULL] = convpolynomialtensor(V,F,O,P)
% [ZEROTH] = convpolynomialtensor(V,F,O,P,ZEROCOMPONENT)
% computes the convolution polynomial F of order O
% of a general Fourier series V in tensor form
% P contains the parameters needed to evaluate the polynomial (if they are needed)
% works in all dimensions; does not assume any symmetry
% the second, optional output FVFULL is the "extended tensor"
% If O is 0 or when F returns a scalar, then the size of FVFULL equals the size of FV.
% If ZEROCOMPONENT is true, only the fourier component with index 0 is returned

sz = size(v);
dim = length(sz);

N = (sz-1)/2;

if exist('parameters','var') && ~isempty(parameters)
      lengthoutput=length(poly(v(:),parameters));
else
      lengthoutput=length(poly(v(:)));
end

s=cell(dim,1);
if order==0 || lengthoutput==1
    % special case of constant polynomial
    fv=0*v;
    for j=1:dim
        s{j}=N(j)+1;
    end
    if exist('parameters','var') && ~isempty(parameters)
        fv(s{1:dim})=poly(0,par);
    else
        % no parameters needed
        fv(s{1:dim})=poly(0);
    end
    fvfull=fv;
    return
end

% lengths of ffts are powers of 2, particularly for intlab
% include zero padding to combat/prevent aliassing
M = 2.^nextpow2(2*order*N+1);  
if exist('zerocomponent','var') && zerocomponent
    % only zeroth component needs to be accurate: half the zero padding needed
    M = 2.^nextpow2(order*N+1);
elseif nargout==1 
    % only components up to original size need to be accurate: 
    % only fraction (order+1)/(2*order) of zero padding needed  
    M = 2.^nextpow2((order+1)*N+1); 
end

v1 = altzeros(M(1:dim),v(1));
M2 = M/2;
M2(M==1) = 0; % take care of dimensions with just one element (N=0)
s1 = M2 + 1 - N;
s2 = M2 + 1 + N;
s1f = M2 + 1 - order*N;
s2f = M2 + 1 + order*N;
sf=cell(dim,1);
fl=cell(dim,1);
for j=1:dim
  s{j} = s1(j):s2(j);
  sf{j} = s1f(j):s2f(j);
  fl{j} = [M2(j)+1:M(j),1:M2(j)];
end

v1(s{1:dim}) = v;
v2 = v1(fl{1:dim}); %fftshift

% fft
w = altfftn(v2);
% apply polynomial
if exist('parameters','var') && ~isempty(parameters)
    wp = reshape(poly(w(:),parameters),size(w));
else
    % no parameters needed
    wp = reshape(poly(w(:)),size(w));
end
if exist('zerocomponent','var') && zerocomponent
    % no ifft required, just the 0-th component    
    fv=sum(wp(:))/prod(M);
else
    vp2 = altifftn(wp);
    vp1 = vp2(fl{1:dim}); %fftshift
    fv = vp1(s{1:dim});
    if nargout==2
        fvfull = vp1(sf{1:dim});
    end
end

return



