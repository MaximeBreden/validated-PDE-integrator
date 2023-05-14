function [polynomials,symmetry]=preprocessnonlinearities(f,epsfactor)
% extract coefficients of polynomials
% and symmetry='even' if nonlinearities respect left-right symmetry
%
% both f=[] and "no input" are fine when there is only
% the leading order derivative and no additional term (e.g. heat equation)
%
% the format of each polynomial is the matlab-standard vector of coefficients 
%
% if interval arithmetic is available then 
% the coefficients of the polynomial are intervals.
% All zero coefficients have zero radius
% but the radii of other intervals take into account rounding of possible nonfloats
% (also for nonzero integer coefficients).
% If you want bigger radii then use epsfactor (default is 1)
% while you may set epsfactor=0 if you know that
% all monomial coefficients are representable as floats (e.g. integers)

global intervalarithmeticavailable 
 
if intervalarithmeticavailable 
   intervalarithmetic=true; 
else
   intervalarithmetic=false; 
end

if ~exist('f','var')
    f = [];
end
    
if ~exist('epsfactor','var') || isempty(epsfactor)
    epsfactor=1;
end

P=length(f);
polynomials=cell(P,1);

% extract coefficients

polynomialtrivial=false(P,1);
for p=1:P
    if ~isempty(f{p})
        % extract coefficients
        polynomials{p}=sym2poly(f{p});
        if intervalarithmetic 
            % find maximum rounding error
            radii=sup(intval(epsfactor)*eps(polynomials{p}));
            % find zero coefficients
            zerocoefficients=(polynomials{p}==0);
            radii(zerocoefficients)=0;
            % turn into intervals with appropriate radii
            polynomials{p}=midrad(polynomials{p},radii);            
        end
    else
        polynomialtrivial(p)=true;
    end
end

% check for possible symmetries

canbeeven=all(polynomialtrivial(2:2:P));
canbeodd=true;
for p=0:P-1
	if ~polynomialtrivial(p+1)
		polyp=fliplr(polynomials{p+1});
		if ~all(polyp(mod(p,2)+1:2:end)==0)
			canbeodd=false;
		end
	end
end
if canbeeven && canbeodd
	symmetry='evenodd';
elseif canbeeven
	symmetry='even';
elseif canbeodd
	symmetry='odd';
else
	symmetry='no';
end

end
