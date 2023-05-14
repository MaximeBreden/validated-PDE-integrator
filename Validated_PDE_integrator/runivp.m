function [x,problem,success,rsol,finaltime,bounds,errorbound] = runivp(filename,directory)
% run the initial value problem with the name filename
% which has initial data in the file ivpdatafilename
% 
% save results (data and figure) in the directory (default 'results')

global intervalarithmeticavailable

if ~exist('directory','var') 
    dirname='results/';
elseif isempty(directory)
    dirname=[];
else
    dirname=[directory,'/'];
end

% get data
datafile=str2func(['ivpdata',filename]);
[problem,plotdata] = datafile();

% floats or intervals or both
if ~isfield(problem,'useintervals') || isempty(problem.useintervals) || ...
                                        ~islogical(problem.useintervals)
    problem.useintervals=[];
    tryboth=true;
else
    tryboth=false;
end
tryintervals = (tryboth || problem.useintervals);
tryfloats = (tryboth || ~problem.useintervals);
if ~intervalarithmeticavailable && ~tryboth && tryintervals
    disp('Interval arithmetic is not available!')
    disp('Continue with floats instead.')
    tryintervals=false;
    tryfloats=true;
end

% solve numerically
[x,problem] = ivpsolve(problem);

% figure
plotsolution(x,problem,plotdata,[dirname,'example',filename]);
save([dirname,'data',filename,'.mat'],'x','problem');

%start proofs
if tryfloats
    % using floats
    useintervals=false;
    tic
    [success,rsol,finaltime,bounds,errorbound] = ivpproof(x,problem,useintervals);
    toc
    tictoc=toc;
    if success 
        save([dirname,'dataprooffloat',filename,'.mat'],'x','problem',...
              'success','rsol','bounds','errorbound','finaltime','useintervals','tictoc');
    else
        disp('Proof with floats failed')
        tryintervals=false;
    end
end

if ~intervalarithmeticavailable && tryintervals && tryboth
    disp('Interval arithmetic is not available!')
    disp('Cannot continue with intervals.')
    tryintervals=false;
end

% now using interval arithmetic
if tryintervals
    useintervals=true;
    tic
    [success,rsol,finaltime,bounds,errorbound] = ivpproof(x,problem,useintervals);
    toc
    tictoc=toc;
    if success 
        save([dirname,'dataproofintval',filename,'.mat'],'x','problem',...
              'success','rsol','bounds','errorbound','finaltime','useintervals','tictoc');
    else
        disp('Proof with intervals failed')
    end
end

end