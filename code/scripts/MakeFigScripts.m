%% sim inference (ie, params are assumed known)

V.datasets  = 13;
V.filters   = [1 4];
V.name      = 'fast_wiener_inf';
V.save      = 1;

P.lam       = 1;
P.sig       = 0.2;

MakeFigs(V,P);

%% sim inference and learning (ie, params are not assumed known)

V.datasets  = 13;
V.filters   = [1.5 4.5];
V.name      = 'fast_wiener_learn';
V.save      = 1;
V.T         = 1000;
V.xlims     = [1 800];
V.dt        = 1/60;
P.lam       = 1;
P.sig       = 0.4;
P.gam       = (1-V.dt/1);
V.fast_iter_max = 3;
V.fast_plot=1;
MakeFigs(V,P);

%% vitro inference and learning (ie, params are not assumed known)

clear, clc
V.datasets  = 4;
V.filters   = [1.5 4.5];
V.name      = 'fast_wiener_vitro';j
V.save      = 1;

MakeFigs(V);

%% vitro burst inference and learning (ie, params are not assumed known)

clear, clc,
V.datasets  = 12;
V.filters   = [1.5 4.5];
V.name      = 'fast_wiener_burst';
V.save      = 1;
V.xlims     = [1100 1700];

MakeFigs(V);

%% smc vs. fast vitro 

clear, clc,
V.datasets  = 12;
V.filters   = 0;
V.name      = 'fast_smc_vitro';
V.save      = 0;
V.fast_do   = 1;
V.smc_do    = 0;
V.preprocess = 0;

MakeFigs(V);
