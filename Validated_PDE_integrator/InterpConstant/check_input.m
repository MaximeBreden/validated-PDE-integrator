function X = check_input(X,dx)

% Keeping only the columns of X compatible with ||f'||<=1. That is, we
% exclude inputs for which |f(x_{k+1}) - f(x_k)| > |x_{k+1}-x_k| for at
% least one k.

dy = abs(X(3:end,:)-X(2:end-1,:));
acceptable_input = not( dy > dx );
acceptable_input = min( acceptable_input, [], 1 ); %combining all the k's
X = X(:,acceptable_input);

end