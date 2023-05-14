function out = obj_func(X,K,M,x)

% The objective function to be minimized.
% P, upper and lower respectively correspond to P_K(y)(t), 
% \overline{\varphi}_K(y)(t) and \underline{\varphi}_K(y)(t) in the paper

t = X(1,:);
y = X(2:end,:);
P = M*y;
upper = -Inf * ones(size(t));
lower = Inf * ones(size(t));
for k = 1 : K/2
    upper = max( upper, min( y(k,:)+t-x(k,:), y(k+1,:)-t+x(k+1,:) ) ); 
    lower = min( lower, max( y(k,:)-t+x(k,:), y(k+1,:)+t-x(k+1,:) ) );
end
if mod(K,2)
    k = (K-1)/2;
    upper = max( upper, min( y(k+1,:)+t-x(k+1,:), -t ) );
    lower = min( lower, max( y(k+1,:)-t+x(k+1,:), +t ) );
end
TK = cos((0:K)' * acos(t));
Pt = sum(P.*TK,1);
out = min( Pt - upper, lower - Pt );

end




