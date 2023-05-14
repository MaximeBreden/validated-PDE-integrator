function normA = matrixnormL1(A,nu1,nu2)
% matrixnormL1(A) computes the matrix normfrom l^1 to l^1
% matrixnormL1(A,nu) computes the matrix norm from l^1_nu to l^1_nu
% matrixnormL1(A,nu1,nu2) computes the matrix norm from l^1_nu1 to l^1_nu2

if ~exist('nu1','var')
    nu1=1;
end
if ~exist('nu2','var')
    nu2=nu1;
end

N1 = (size(A,2)-1)/2; 
N2 = (size(A,1)-1)/2; 
n1 = (-N1:N1); 
n2 = (-N2:N2); 
normA = max(((nu2.^abs(n2))*abs(A))./(nu1.^abs(n1)));

end
