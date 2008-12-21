% this function takes in either no input (in which case it simulates data),
% or 1 input (which is the fluorescence time series. in either case, it
% runs some OOPSI filters to infer stuff

% set simulation metadata
Sim.T   = 100;
Sim.dt      = 0.031;                                % median frame duration
Sim.MaxIter  = 8;                                  % max # of EM iterartions
Sim.Plot     = 1;                                   % plot results with each iteration

% initialize FOOPSI parameters
P.tau   = 2;                                        % calcium decay time constant (sec)
P.sig   = 0.1;                                      % standard deviation of noise (\mu M)
Nsp     = 50;                                       % expected spike rate
P.lam   = Sim.T/(Nsp*Sim.dt);                       % expected jump size ber time bin


% if no input fluorescence vector    
    n=rand(Sim.T,1)<P.lam*Sim.dt*exp(-P.lam*Sim.dt);
    C=filter(1,[1 -(1-Sim.dt/P.tau)],n);
    F=C+P.sig.*randn(Sim.T,1);

% infer spikes and estimate parameters
Algs = [14];
for m=Algs;
    Sim.Alg = m;
%     P.lam   = 1;
P.sig=P.sig*2;
    I{m}    = DataComp14(F,P,Sim);
end

%% plot results
figure(1); clf,
row = 3;
nAlgs = numel(Algs);
nrows = row+nAlgs; 
subplot(nrows,1,1), plot(F), axis('tight'), ylabel('F')
if nargin==0
    subplot(nrows,1,2), plot(C), axis('tight'), ylabel('C')
    subplot(nrows,1,3), bar(n), axis('tight'), ylabel('n')
end
for m=Algs
    row=row+1;
    subplot(nrows,1,row), bar(I{m}.n), axis('tight'), ylabel(I{m}.name)
end