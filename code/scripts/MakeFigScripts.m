%% make figure 2

V.datasets  = 13;
V.filters   = [1 4];
V.name      = 'fast_wiener_inf';
V.save      = 1;

P.lam       = 1;
P.sig       = 0.2;

MakeFigs(V,P);

%% make figure 4

V.datasets  = 13;
V.filters   = [1 4];
V.name      = 'fast_wiener_learn';
V.save      = 1;

P.lam       = 1;
P.sig       = 0.2;

MakeFigs(V,P);
