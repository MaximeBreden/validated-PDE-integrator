function maxdegree = nonlinmaxdegree(problem)
% maximum degree of all the nonlinearities

maxdegree=0;
for p=1:length(problem.pde.polynomials)
    maxdegree=max(maxdegree,nonlindegree(p-1,0,problem));
end

end
