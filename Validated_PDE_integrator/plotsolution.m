function plotsolution(x,problem,plotdata,filename)
% plots the (numerical) solution x of the initial value problem
% where x is of size (2N+1)x(K+1)xD
%
% if filename is provided and nonempty then the plot is saved
%
% plotdata may contain
% plotdata.viewangle : then view(viewangle) is used
% plotdata.numpoints : then (dafault is [101,101])
% the number of gridpoints in time is numpoints(1) and
% the number of gridpoints in space is numpoints(2)

if ~exist('plotdata','var') 
    plotdata=[];
end

if ~isfield(plotdata,'viewangle') || isempty(plotdata.viewangle)
    plotdata.viewangle=3; % default in 3D
end
if ~isfield(plotdata,'numpoints') || isempty(plotdata.numpoints)
    plotdata.numpoints=[101,101];
end

% number of points used in the plot
timepoints=plotdata.numpoints(1);
spacepoints=plotdata.numpoints(2);

N=(size(x,1)-1)/2;
D=size(x,3);
tgrid=problem.timegridfloats;

% adjust number of points in time to align with a possibly uniform grid
%timepoints=D*ceil(timepoints/D)+1;
times=linspace(0,tgrid(end),timepoints);
space=linspace(0,2*pi,spacepoints)';
solution=zeros(spacepoints,timepoints);

for d=1:D
    % time points beloning to d-th domain
    td=(times>tgrid(d) & times<=tgrid(d+1));
    if d==1
        % include left endpoint for first domain
        td(1)=true;
    end
    % rescale time to [-1,1]
    s=2*((times(td)-tgrid(d))/(tgrid(d+1)-tgrid(d)))-1;
    % evaluate in time
    sol=cheb_eval(x(:,:,d),s);
    % evaluate in space
    sol=exp(1i*space*(-N:N))*sol;
    % remove spurious imaginary parts
    sol=real(sol);
    % add the solution on domain d to the full solution
    solution(:,td)=sol;
end

% plot
maxtime=problem.scalefloats.time*tgrid(end);
maxspace=problem.scalefloats.space*2*pi;
plottimes=linspace(0,maxtime,timepoints);
plotspace=linspace(0,maxspace,spacepoints);
[xvar,tvar]=ndgrid(plotspace,plottimes);

fig=figure;
surf(xvar,tvar,solution);

% some print setting for pdf
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

xlabel('$x$','Interpreter','latex','FontSize',24)
ylabel('$t$','Interpreter','latex','FontSize',24)
zlabel('$u$','Interpreter','latex','FontSize',24)
axis tight
view(plotdata.viewangle)
drawnow

% save plot
if exist('filename','var') && ~isempty(filename)
    print(fig,filename,'-dpdf','-r300')
end

end

