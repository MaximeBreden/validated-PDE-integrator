% This produces all figures (the numbering of the figures corresponds to
% the one used in the paper)

% directory to get the data from
datadir='data'; 
% directory to put the figures in
figuredir='plots';
mkdir(figuredir);

%%%%%%%%%%%%
% Figure 1 %
%%%%%%%%%%%% 
problemname='SwiftHohenberg1';
datafilename=[datadir,'/dataintval',problemname];
load(datafilename,'x','problem');
figurefilename=[figuredir,'/example',problemname];
plotsolution(x,problem,[],figurefilename);

plotdata.viewangle = [-160 30];
figurefilename=[figuredir,'/example',problemname,'a'];
plotsolution(x,problem,plotdata,figurefilename);
plotdata.viewangle = [];

%%%%%%%%%%%%
% Figure 2 %
%%%%%%%%%%%%
problemname='Fisher1';
datafilename=[datadir,'/dataintval',problemname];
load(datafilename,'x','problem');
figurefilename=[figuredir,'/example',problemname];
plotsolution(x,problem,[],figurefilename);

%%%%%%%%%%%%
% Figure 3 %
%%%%%%%%%%%%
problemname='SwiftHohenberg2';
datafilename=[datadir,'/dataintval',problemname];
load(datafilename,'x','problem');
figurefilename=[figuredir,'/example',problemname];
plotsolution(x,problem,[],figurefilename);

%%%%%%%%%%%%
% Figure 4 %
%%%%%%%%%%%%
problemname='Kuramoto1';
datafilename=[datadir,'/dataintval',problemname];
load(datafilename,'x','problem');
figurefilename=[figuredir,'/example',problemname];
plotsolution(x,problem,[],figurefilename);

%%%%%%%%%%%%
% Figure 5 %
%%%%%%%%%%%%
problemname='OhtaKawasaki';
datafilename=[datadir,'/dataintval',problemname];
load(datafilename,'x','problem');
figurefilename=[figuredir,'/example',problemname];
plotsolution(x,problem,[],figurefilename);

%%%%%%%%%%%%
% Figure 6 %
%%%%%%%%%%%%
problemname='Fisher4';
datafilename=[datadir,'/timestepproofintval',problemname];
load(datafilename,'problem','floats');
figurefilename=[figuredir,'/example',problemname];
plotdata.scale=problem.scalefloats;
intermediate=floats.intermediate;
times=floats.times;
plottimesteps(times,intermediate,plotdata,figurefilename)
