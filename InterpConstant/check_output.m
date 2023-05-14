function [X, fX, current_min] = check_output(X,fX,threshold)

% We keep only the columns of X for which f could be below threshold, and
% sort those columns Xj by increasing value of the inf of f(Xj). Hence, if
% we later decrease the thresold and get thresold < f(Xj) for some j, we 
% know we can discard all the columns having indices > j.

current_min = min(fX);
threshold = min(threshold,sup(current_min));
acceptable_output = not( inf(fX) > threshold ); 
X = X(:,acceptable_output);
fX = fX(acceptable_output);
[~,I] = sort(inf(fX));
X = X(:,I);
fX = fX(I);

end