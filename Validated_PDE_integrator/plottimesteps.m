function plottimesteps(times,intermediate,plotdata,filename)
% plots the intermediate profiles of a time stepping solution
% both as 3D plots (in two ways) and as a 2D plot
%
% if filename is provided and nonempty then the plots are saved
%
% plotdata may contain
% plotdata.scale.space (default 1) 
% plotdata.scale.time (default 1) 
% plotdata.viewangle : then view(viewangle) is used
% plotdata.numpoints (default is [101,101]) :
% the number of gridpoints in time is numpoints(1) and
% the number of gridpoints in space is numpoints(2)

if ~exist('plotdata','var') 
    plotdata=[];
end

if ~isfield(plotdata,'scale')
    plotdata.scalefloat.space=1;
    plotdata.scalefloat.time=1;
end
if ~isfield(plotdata,'viewangle') || isempty(plotdata.viewangle)
    plotdata.viewangle=3; % default in 3D
end
if ~isfield(plotdata,'numpoints') || isempty(plotdata.numpoints)
    plotdata.numpoints=[101,101];
end

% number of grid points used in the plot
spacepoints=plotdata.numpoints(2);

fourier=altmid(intermediate.data);
times=altmid(times);
N=(size(fourier,1)-1)/2;
M=length(times);

space=linspace(0,2*pi,spacepoints)';
profiles=real(exp(1i*space*(-N:N))*fourier);

plotspace=linspace(0,plotdata.scale.space*2*pi,spacepoints);
timescale=plotdata.scale.time;

fig3dslice=figure;
for m=1:M
    plot3(plotspace,repmat(timescale*times(m),spacepoints,1),...
                                profiles(:,m),'b','LineWidth',2);
    hold on
end
fig2d=figure;
for m=1:M
    plot(plotspace,profiles(:,m));
    hold on
end

fig3d=figure;
[xvar,tvar]=ndgrid(plotspace,timescale*times);
surf(xvar,tvar,profiles);

% some print setting for pdf
set(fig3d,'Units','Inches');
pos3d = get(fig3d,'Position');
set(fig3d,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos3d(3), pos3d(4)])
set(fig3dslice,'Units','Inches');
pos3dslice = get(fig3dslice,'Position');
set(fig3dslice,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos3dslice(3), pos3dslice(4)])
set(fig2d,'Units','Inches');
pos2d = get(fig2d,'Position');
set(fig2d,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos2d(3), pos2d(4)])

% finish 3D figure
figure(fig3d)
hold off
xlabel('$x$','Interpreter','latex','FontSize',24)
ylabel('$t$','Interpreter','latex','FontSize',24)
zlabel('$u$','Interpreter','latex','FontSize',24)
axis tight
view(plotdata.viewangle)
drawnow

% finish 3D slice figure
figure(fig3dslice)
hold off
xlabel('$x$','Interpreter','latex','FontSize',24)
ylabel('$t$','Interpreter','latex','FontSize',24)
zlabel('$u$','Interpreter','latex','FontSize',24)
axis tight
view(plotdata.viewangle)
drawnow

% finish 2D figure
figure(fig2d)
hold off
xlabel('$x$','Interpreter','latex','FontSize',20)
ylabel('$u$','Interpreter','latex','FontSize',20)
axis tight
drawnow

% save plots 
if exist('filename','var') && ~isempty(filename)
    print(fig3d,[filename,'-3D.pdf'],'-dpdf','-r300')
    print(fig3dslice,[filename,'-3Dslice.pdf'],'-dpdf','-r300')
    print(fig2d,[filename,'-2D.pdf'],'-dpdf','-r300')
end

end

