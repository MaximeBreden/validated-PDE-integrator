% This file runs all proofs of the theorems

initialize

% change to true only when running on a big machine (we tested with 256Gb of RAM)
lotsofmemoryavailable = false;

% directory to get the data from
datadir='data'; 
% directory to put the results
resultdir='reproduce';
mkdir(resultdir);

if ~intervalarithmeticavailable
    % run all "proofs" with floats
    useintervals=false;
    datatype='float';
else
    % run all proofs rigorously
    useintervals=true;
    datatype='intval';
end

%%%%%%%%%%%%%%%
% Theorem 1.5 %
%%%%%%%%%%%%%%%

if lotsofmemoryavailable
    problemname='SwiftHohenberg1';
    datafilename=[datadir,'/data',datatype,problemname];
    load(datafilename,'x','problem');
    fprintf("\n\nRunning the 'accurate' proof for the Swift-Hohenberg equation (Theorem 1.5)\n\n") 
    tic
    [success,rsol,finaltime,bounds,errorbound] = ivpproof(x,problem,useintervals);
    toc
    tictoc=toc;
    resultsfilename=[resultdir,'/dataproof',datatype,problemname];
    save(resultsfilename,'x','problem','success','rsol','bounds','errorbound','finaltime','useintervals','tictoc');
end

%%%%%%%%%%%%%%%
% Theorem 7.1 %
%%%%%%%%%%%%%%%

problemname='Fisher1';
datafilename=[datadir,'/data',datatype,problemname];
load(datafilename,'x','problem');
fprintf("\n\nRunning the proof for Fisher's equation (Theorem 7.1)\n\n")
tic
[success,rsol,finaltime,bounds,errorbound] = ivpproof(x,problem,useintervals);
toc
tictoc=toc;
resultsfilename=[resultdir,'/dataproof',datatype,problemname];
save(resultsfilename,'x','problem','success','rsol','bounds','errorbound','finaltime','useintervals','tictoc');

%%%%%%%%%%%%%%%
% Theorem 7.3 %
%%%%%%%%%%%%%%%

problemname='SwiftHohenberg2';
datafilename=[datadir,'/data',datatype,problemname];
load(datafilename,'x','problem');
fprintf("\n\nRunning the 'easy' proof for the Swift-Hohenberg equation (Theorem 7.3)\n\n")
tic
[success,rsol,finaltime,bounds,errorbound] = ivpproof(x,problem,useintervals);
toc
tictoc=toc;
resultsfilename=[resultdir,'/dataproof',datatype,problemname];
save(resultsfilename,'x','problem','success','rsol','bounds','errorbound','finaltime','useintervals','tictoc');

%%%%%%%%%%%%%%%
% Theorem 7.5 %
%%%%%%%%%%%%%%%

if lotsofmemoryavailable
    problemname='Kuramoto1';
    datafilename=[datadir,'/data',datatype,problemname];
    load(datafilename,'x','problem');
    fprintf("\n\nRunning the proof for the Kuramoto-Sivashinsky equation (Theorem 7.5)\n\n")
    tic
    [success,rsol,finaltime,bounds,errorbound] = ivpproof(x,problem,useintervals);
    toc
    tictoc=toc;
    resultsfilename=[resultdir,'/dataproof',datatype,problemname];
    save(resultsfilename,'x','problem','success','rsol','bounds','errorbound','finaltime','useintervals','tictoc');
end

%%%%%%%%%%%%%%%
% Theorem 7.6 %
%%%%%%%%%%%%%%%

if lotsofmemoryavailable
    problemname='OhtaKawasaki';
    datafilename=[datadir,'/data',datatype,problemname];
    load(datafilename,'x','problem');
    fprintf("\n\nRunning the proof for the Ohta-Kawasaki equation (Theorem 7.5)\n\n")
    tic
    [success,rsol,finaltime,bounds,errorbound] = ivpproof(x,problem,useintervals);
    toc
    tictoc=toc;
    resultsfilename=[resultdir,'/dataproof',datatype,problemname];
    save(resultsfilename,'x','problem','success','rsol','bounds','errorbound','finaltime','useintervals','tictoc');
end

%%%%%%%%%%%%%%%
% Theorem 7.9 %
%%%%%%%%%%%%%%%

if intervalarithmeticavailable
    problemname='Fisher4';
    datafilename=[datadir,'/timestepproof',datatype,problemname];
    load(datafilename,'problem');
    fprintf("\n\nRunning the proof for the Fisher'e equation with time-stepping (Theorem 7.9)\n\n")
    tic
    [final,times,intermediate,uniformerrorbound] = timestepping(problem,useintervals);
    toc
    tictoc=toc;
    resultsfilename=[resultdir,'/timestepproof',datatype,problemname];
    save(resultsfilename,'problem','final','times','intermediate','uniformerrorbound','useintervals','tictoc');
else
    disp('One cannot properly mimic the time stepping proof without interval arithmetic.')
end
