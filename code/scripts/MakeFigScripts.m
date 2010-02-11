%% make figure 2

V.datasets  = 13;
V.filters   = [1 4];
V.name      = 'fast_wiener_sim';
V.save      = 1;

P.lam       = 1;
P.sig       = 0.2;

MakeFigs(V,P);