function [final,problem,success,times,intermediate,uniformerrorbound] = runtimestepping(filename,directory)
% run the time stepping problem with the name filename
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

%start time stepping
if tryfloats
    % using floats
    useintervals=false;
    tic
    [final,times,intermediate,uniformerrorbound] = timestepping(problem,useintervals);
    toc
    tictoc=toc;
    if ~any(isnan(final.data))
        success=true;
        floats.intermediate=intermediate;
        floats.times=times;
        save([dirname,'timestepprooffloat',filename,'.mat'],'problem','success',...
            'final','times','intermediate','uniformerrorbound','useintervals','tictoc','floats');
        disp('Time stepping with floats was successful')
        disp(['The uniform error bound is ',num2str(uniformerrorbound)])
        disp('But this is not mathematically rigorous')
    else
        success=false;
        disp('Time stepping with floats failed')
        disp(['Could verify ',num2str(length(times)-1),' steps out of ',...
                                    num2str(length(problem.timegrid)-1)]);
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
    [final,times,intermediate,uniformerrorbound] = timestepping(problem,useintervals);
    toc
    tictoc=toc;
    if ~any(isnan(final.data))
        success=true;
        % save floats for plotting purposes
        floats.intermediate=intermediate;
        floats.intermediate.data=altmid(floats.intermediate.data);
        floats.times=altmid(times);
        save([dirname,'timestepproofintval',filename,'.mat'],'problem','success',...
            'final','times','intermediate','uniformerrorbound','useintervals','tictoc','floats');
        disp('Time stepping with intervals was successful')
        disp(['The uniform error bound is ',num2str(uniformerrorbound)])
        disp('Full proof including interval arithmetic')
    else
        success=false;
        disp('Time stepping with intervals failed')
        disp(['Could verify ',num2str(length(times)-1),' steps out of ',...
                                    num2str(length(problem.timegrid)-1)]);
    end
end

% plot
plotdata.scale=problem.scalefloats;
plottimesteps(floats.times,floats.intermediate,plotdata,[dirname,'timesteps',filename]);

end