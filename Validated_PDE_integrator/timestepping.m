function [final,times,intermediate,uniformerrorbound] = timestepping(problem,useintervals)
% performs time stepping for problem until it fails
%
% problem contains all the information of the IVP
% as well as the data for solving it 
% (timegrid, proof constants, etc.)
% 
% final contains the data and tailbound for the final 
% if successful until the end (otherwise NaN)
% times contains the successful times (a subset of problem.timegrid)
% and intermediate.data (and intermediate.tailbounds) at these times
%
% if useintervals is true then all computations will be rigorous
% (unless intervalarithmetic is not available)

%%%%%%%%%%%%%%%%%%
% Initialization %
%%%%%%%%%%%%%%%%%%

global intervalarithmeticavailable
if intervalarithmeticavailable && useintervals
    disp('Including interval arithmetic')
    timevector=intval(problem.timegrid);
else
   if ~intervalarithmeticavailable
        disp('No interval arithmetic available!!')
        disp('Continuing without intervals')
   else
       disp('No interval arithmetic!')
   end
   timevector=altmid(problem.timegrid);
   useintervals=false;
end    

% initialize intermediate points
nsteps=length(problem.timegrid)-1;
N=(length(problem.initial.data)-1)/2;
intermediate.data=altzeros([2*N+1,nsteps+1],useintervals);
intermediate.tailbound=zeros([1,nsteps+1]);
if useintervals
    intermediate.data(:,1)=problem.initial.data;
else
    intermediate.data(:,1)=altmid(problem.initial.data);
end
intermediate.tailbound(1)=problem.initial.tailbound;

% set up data for time stepping
D0=length(problem.subgrid)-1;
step=1;
subproblem=problem;
deltatimes=timevector(2:end)-timevector(1:end-1);
% we continue timestepping as long as we have not failed
failed=false;

%%%%%%%%%%%%%%%%%
% Time stepping %
%%%%%%%%%%%%%%%%%

uniformerrorbound=0;

while ~failed && step<=nsteps
    disp(['time step ' ,int2str(step)]); 
    
    % set up domain of integration
    subproblem.timegrid=problem.subgrid*deltatimes(step);
    subproblem.domains=subproblem.timegrid(2:D0+1)-subproblem.timegrid(1:D0);
    subproblem.proof.rstar=repmat(problem.proof.rstar(step),D0,1);
    subproblem.proof.eps0=repmat(problem.proof.eps0(step),D0,1);
    subproblem.proof.eps1=repmat(problem.proof.eps1(step),D0,1);

    % solve numerically
    [x,subproblem] = ivpsolve(subproblem);
    [success,rsol,finaltime,~,errorbound] = ivpproof(x,subproblem,useintervals);

    if success
        % enclose final time
        [finaldata,finaltailbound] = finaltimeenclosure(x,subproblem,rsol,finaltime,useintervals);
        uniformerrorbound=max(uniformerrorbound,errorbound);

        % and update initial data
        subproblem.initial.data = finaldata;
        subproblem.initial.tailbound = finaltailbound;

        % store and display
        intermediate.data(:,step+1)=subproblem.initial.data;
        intermediate.tailbound(step+1)=subproblem.initial.tailbound;
        disp(['In step ',num2str(step),' the tailbound is ',num2str(subproblem.initial.tailbound)]);
        step=step+1; 
    else
        disp(['Time stepping failed at step ',num2str(step)]);
        failed=true;
    end
end

if failed
    final.data=NaN*problem.initial.data;
    final.tailbound=NaN;
    % truncated data
    times=timevector(1:step);
    intermediate.data=intermediate.data(:,1:step);
    intermediate.tailbound=intermediate.tailbound(1:step);
else
    % successful until the end
    times=timevector;
    final.data=intermediate.data(:,nsteps+1);
    final.tailbound=intermediate.tailbound(nsteps+1);
end

end

