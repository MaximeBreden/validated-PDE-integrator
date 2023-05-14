function [finaldata,finaltailbound] = finaltimeenclosure(x,problem,rsol,finaltime,useintervals)
% find enclosure of solution at final time 
% for time stepping purposes
%
% recursively if D>1
% 
% the input finaltime contains bounds on some integral terms 
% to enclose the solution at the end of the time interval 

N=(size(x,1)-1)/2;
D=size(x,3);
finaltailbound=problem.initial.tailbound;
finaldata=problem.initial.data;

if useintervals        
    nu = intval(problem.proof.nu); 
    weights=nu.^(-abs(-N:N)');
    problem.domains=problem.timegrid(2:D+1)-problem.timegrid(1:D);
    tau=problem.domains/2;
    for d=1:D
        % "recursively" over each subdomain
        rdball=cintval(0,rsol(d));
        rd=intval(rsol(d));

        % initial value term
        lambda=problem.semigroup.Lambda(:,d);
        Q=problem.semigroup.Q(:,:,d);
        Qinv=problem.semigroup.Qinv(:,:,d);
        expL=Q*diag(exp(2*tau(d)*lambda))*Qinv;
        % the finaltimedata from previous interval are initial data
        finaldata=expL*finaldata;

        % integral term
        % careful with forcing interval arithmetic here
        integraldataterm=finaltime.Y.data(:,d) + ...
                            rdball*finaltime.Z.data(:,d) + ...
                            rdball*rdball/2*finaltime.W.data(:,d);
        finaldata=finaldata+integraldataterm;

        % intersect with the uniform bound
        theta=max(max(problem.proof.eps0(d),problem.proof.eps1(d)),1);
        uniformdatabound=2*sum(x(:,:,d),2)-x(:,1,d)+theta*rdball*weights;       
        finaldata=intersect(finaldata,uniformdatabound);
        if isfield(problem,'symmetry') && strcmp(problem.symmetry,'cosineseries') 
            finaldata=real(finaldata);
        end

        % tail term
        % initial data term
        % the finaltimedata from previous interval are initial data
        finaltailbound=boundat1tailexp(N,d,problem)*finaltailbound;
        % integral term
        integraltailterm=finaltime.Y.tailbound(d) + ...
                            rd*finaltime.Z.tailbound(d) + ...
                            rd^2/2*finaltime.W.tailbound(d);
        finaltailbound=finaltailbound+integraltailterm;
        % intersect with uniform bound
        uniformtailbound=rd*problem.proof.eps1(d);
        finaltailbound=min(finaltailbound,uniformtailbound);
        finaltailbound=altsup(finaltailbound);
    end
else
    % no internal arithmetic
    finaltimedata1=finaltime.Y.boundary+finaltime.Y.data(:,D);
    finaltimedata2=2*sum(x(:,:,D),2)-x(:,1,D);
    % two ways to get the approximate data
    finaldata=(finaltimedata1+finaltimedata2)/2;

    for d=1:D
        finaltailbound=altsup(boundat1tailexp(N,d,problem))*finaltailbound;
        integraltailterm=finaltime.Y.tailbound(d) + ...
                            rsol(d)*finaltime.Z.tailbound(d) + ...
                            rsol(d)^2/2*finaltime.W.tailbound(d);
        uniformtailbound=rsol(d)*problem.proof.eps1(d);
        finaltailbound=min(finaltailbound,uniformtailbound);
        finaltailbound=finaltailbound+integraltailterm;
    end
end

end