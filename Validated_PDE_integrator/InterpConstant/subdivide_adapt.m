function Y = subdivide_adapt(X,K,nb_sub,lmax)

% Given X = [X1,...,X_N2], which we view as N2 vectors of intervals, the
% output Y = [Y1,...,Y_M] is another collection of vectors of intervals, 
% with M>N2, such that the union of the Xi is contained in the union of
% the Yi. The Yi are obtained by subdividing the Xi. The subdivisions are 
% controlled by the parameter nb_sub: each component of each Xi is split
% into nb_sub components (of equal length). The optional input parameter
% lmax can be used to control the maximal number M of Yi produced,
% replacing the inputed nb_sub by a smaller value if necessary.

if mod(K,2)==0
    %if K is even, 0 is an interpolation node at which the value of f is
    %fixed to 0, so there is no need to subdivide [0,0].
    X = X([1:K/2+1,K/2+3:end],:);
end

[N1,N2] = size(X);

if nargin == 4
    nb_sub = max(2,min(nb_sub,floor((lmax/N2)^(1/N1))));
%     if nb_sub == 1
%         disp("no further subdivision here")
%     end
end
Y = intval(zeros(N1,N2*(nb_sub^N1)));
resc = linspace(0,1,nb_sub+1)';
for n = 1:N1
    Xn = inf(X(n,:)) + resc*diam(X(n,:)); % subdivision of the n-th components
    Xn(1,:) = inf(X(n,:));
    Xn(end,:) = sup(X(n,:)); % making sure that we did not loose small bits at the extremities due to floatting point errors
    Xn = infsup(Xn(1:end-1,:),Xn(2:end,:));
    Xn = repmat(repelem(Xn,nb_sub^(n-1),1),nb_sub^(N1-n),1); % replicating the subintervals according to the number of components
    Xn = reshape(Xn, [1,numels(Xn)]);
    Y(n,:) = Xn;
end

if mod(K,2)==0
    %putting back the zeros
    Y = [Y(1:K/2+1,:);zeros(1,size(Y,2));Y(K/2+2:end,:)];
end