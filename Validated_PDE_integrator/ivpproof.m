function [success,rsol,finaltime,bounds,errorbound]=ivpproof(x,problem,intervals)
% proof with or without intervalarithmetic
% if x(1) is an interval or intervals==true 
% then interval arithmetic is used if it is available

global intervalarithmeticavailable 

% decide whether to use intervals or not
useintervals=false;
if exist('intervals','var')
    if intervals==true 
        if intervalarithmeticavailable
            useintervals=true;
        else
            disp(['No interval arithmetic available!! '...
                        'Continuing without intervals'])
        end

    end
elseif intervalarithmeticavailable && altisintval(x(1))
    % default is to use intervals if x(1) is an interval
    useintervals=true;
end

if useintervals
    x=intval(x);
    disp('Full proof including interval arithmetic')
else
    x=altmid(x);
    disp('Proof attempt without interval arithmetic')
end

if ~isfield(problem,'semigroup') 
    problem.semigroup.type='elaborate';
end
% interval arithmetic semigroup data if x are intervals
problem.semigroup=fixsemigroup(x,problem);

% run the proof
[success,rsol,finaltime,bounds,errorbound]=provesol(x,problem);

end