% This produces all data, theorems and figures

initialize
resultdir='reproduce';
mkdir(resultdir);

% change to true only when running on a big machine (we tested with 256Gb of RAM)
lotsofmemoryavailable = false;

%%%%%%%%%%%%%%
% Fisher-KPP %
%%%%%%%%%%%%%%

fprintf("\nGenerating data and running the proof for Fisher's equation (Theorem 7.1)\n\n")
% This computes the data and runs the proof of Theorem 7.1 and Figure 2
runivp('Fisher1',resultdir);

fprintf("\n\nGenerating data and running the proof for Fisher's equation with slightly more subdomains (Remark 7.2)\n\n")
% This runs a similar computation but with more accuracy (Remark 7.2)
runivp('Fisher2',resultdir);

fprintf("\n\nGenerating data and running the proof for Fisher's equation with a naive choice of L (Section 7.2)\n\n")
% This runs the proof with a naive choice for the semigroup (the operator L)
% see Section 7.2 of the paper
runivp('Fisher3',resultdir);

fprintf("\n\nGenerating data and running the proof for Fisher's equation with time stepping (Theorem 7.9)\n\n")
% This runs the time stepping problem (rather than domain decomposition)
% see Section 7.3, Theorem 7.9 and Figure 6
runtimestepping('Fisher4',resultdir);

%%%%%%%%%%%%%%%%%%%
% Swift-Hohenberg %
%%%%%%%%%%%%%%%%%%%

if lotsofmemoryavailable 
    fprintf("\n\nGenerating data and running the 'accurate' proof for the Swift-Hohenberg equation (Theorem 1.5)\n\n")
    % This computes the data and runs the proof of Theorem 1.5 and Figure 1
    runivp('SwiftHohenberg1',resultdir);
end

fprintf("\n\nGenerating data and running the 'easy' proof for the Swift-Hohenberg equation (Theorem 7.3)\n\n")
% This computes the data and runs the proof of Theorem 7.3 and Figure 3
runivp('SwiftHohenberg2',resultdir);

%%%%%%%%%%%%%%%%%%%%%%%%
% Kuramoto-Sivashinsky %
%%%%%%%%%%%%%%%%%%%%%%%%

if lotsofmemoryavailable 
    fprintf("\n\nGenerating data and running the proof for the Kuramoto-Sivashinsky equation (Theorem 7.5)\n\n")
    % This computes the data and runs the proof of Theorem 7.5 and Figure 4
    runivp('Kuramoto1',resultdir);
end

fprintf("\n\nGenerating data and running an 'easier' proof for the Kuramoto-Sivashinsky equation (not in the paper, same solution as in Theorem 7.5, but with a worse error bound)\n\n")
% This computes the data and runs the proof with less accuracy
% but still producing a proof
runivp('Kuramoto2',resultdir);

%%%%%%%%%%%%%%%%%
% Ohta-Kawasaki %
%%%%%%%%%%%%%%%%%

if lotsofmemoryavailable 
    fprintf("\n\nGenerating data and running the proof for the Ohta-Kawasaki equation (Theorem 7.6)\n\n")
    % This computes the data and runs the proof of Theorem 7.6 and Figure 5
    runivp('OhtaKawasaki',resultdir);
end


