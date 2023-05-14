function cheb = fixcheb(K1,intvaltest)
% stores some Chebyshev data
% with KK1=2.^nextpow2(K1) gridpoints

cheb.K1=K1;
cheb.KK1=2.^nextpow2(K1);
cheb.weights=cheb_quadrature_weights(cheb.KK1,intvaltest);
cheb.grid=cheb_grid(cheb.KK1,intvaltest);

end

