function problem = preprocessdata(problem,initialdata,griddata,symmetry,scale,subgriddata)
% preprocessing of the input data to define
% problem.initial.data
% problem.initial.tailbound
% problem.timegrid
% problem.timegridfloats
% problem.symmetry
% problem.scale
% problem.scalefloats
% problem.subgrid (for time stepping)
%
% keeps track of rounding errors when converting to intervals

global intervalarithmeticavailable

if all(imag(initialdata)==0)
    initialdatareal=true;
else
    initialdatareal=false;
end
if all(real(initialdata)==0)
    initialdataimag=true;
else
    initialdataimag=false;
end

% turn into intervals if appropriate
if intervalarithmeticavailable && ~isintval(initialdata(1))
    % factor sqrt(2) for complex-valued data
    if initialdatareal || initialdataimag
        initialdataradii=eps(initialdata);
    else
        initialdataradii=sup(intval('sqrt2')*eps(initialdata));
    end
    zeroinitialdata=(initialdata==0);
    initialdataradii(zeroinitialdata)=0;
    initialdata=midrad(initialdata,initialdataradii);
end

% symmetrize to make initial data real; note the factor 1/2
initialdata=(initialdata+conj(initialdata(end:-1:1)))/2;

% store initial data
problem.initial.data=initialdata;
problem.initial.tailbound=0; % in the nu-norm

% time grid
D=griddata.D;
T=altmid(griddata.T);
if ~isfield(griddata,'skew')
    griddata.skew=0;
end
if ~isfield(griddata,'dip')
    griddata.dip=0;
end
s=linspace(0,1,D+1);
% define the grid
problem.timegridfloats=T*s.*(1+griddata.skew*(s-1)+griddata.dip*(s-1).*(s-0.5));
if altisintval(griddata.T)
    problem.timegrid=intval(problem.timegridfloats);
    % replace final gridpoint by the integration time as an interval
    problem.timegrid(D+1)=griddata.T;
else
    problem.timegrid=problem.timegridfloats;
end
if ~all((problem.timegrid(2:D+1)-problem.timegrid(1:D))>0) 
    disp('Negative step sizes: need to adjust skew-dip');
    error('Negative step sizes: need to adjust skew-dip');
end

% determine symmetry
if initialdatareal && contains(symmetry,'even')
    problem.symmetry='cosineseries';
elseif initialdataimag && contains(symmetry,'odd')
    problem.symmetry='sineseries';
else
    problem.symmetry='no';
end

% set scales
if ~exist('scale','var') || isempty(scale)
    scale.space=1;
    scale.time=1;
end
problem.scale=scale;
problem.scalefloats.space=altmid(problem.scale.space);
problem.scalefloats.time=altmid(problem.scale.time);

% determine subgrid
if ~exist('subgriddata','var')
    subgriddata.D0=1;
    subgriddata.skew=0;
    subgriddata.dip=0;
end
if ~isfield(subgriddata,'skew')
    subgriddata.skew=0;
end
if ~isfield(subgriddata,'dip')
    subgriddata.dip=0;
end
D0=subgriddata.D0;
s0=linspace(0,1,D0+1);
problem.subgrid=s0.*(1+subgriddata.skew*(s0-1)+subgriddata.dip*(s0-1).*(s0-0.5));
if ~all((problem.subgrid(2:D0+1)-problem.subgrid(1:D0))>0) 
    disp('Negative step sizes: need to adjust skew-dip for subgrid');
    error('Negative step sizes: need to adjust skew-dip for subgrid');
end

end