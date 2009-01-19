%% this file compares various filters for real data.  specifically D080218a
%% trial 24

clear,
foldname    = 'TS108_6';
i=2;
if i<10                                         % get tif file name
    savname=[foldname '_00' num2str(i) '.mat'];
else
    savname=[foldname foldname '_0' num2str(i) '.mat'];
end

load(['/Users/joshyv/Research/projects/oopsi/spatial-filter/' savname])

F=D{2}.F';

% 1) set simulation metadata
Sim.T       = size(F,1);                        % # of time steps
Sim.dt      = 30/Sim.T;                         % time step size
Sim.Plot    = 1;                                % whether to plot
Sim.MaxIter = 0;                                % max number of iterations

siz         = size(F);                          % size of data mat
O           = 1+0*F(1,:);                       % vector of size of image

% 2) initialize parameters
P.a     = .1*O;
P.b     = .2*O;
P.sig   = .03;                                  % stan dev of noise
tau     = 2;                                 % decay time constant
P.gam   = 1-Sim.dt/tau;
P.lam   = 30;                                   % rate-ish, ie, lam*dt=# spikes per second
P.l     = 1e99;                                 % initialize likelihood

% 3) infer spikes and estimate parameters

[I{1}.n I{1}.P] = FOOPSI2_54(F,P,Sim);

%%
% 4) plot results
fig     = figure(1); clf,
nrows   = 3;                                    % set number of rows
h       = zeros(nrows);
Pl.xlims= [5 Sim.T];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
n       = 0*F(:,1);
n(R.spt)= 1;
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
i       = 0;

% plot fluorescence data
i=1; h(1) = subplot(nrows,1,i);
Pl.label = 'Fluorescence';
Pl.color = 'k';
Plot_X(Pl,F);

% plot spike train
i=i+1; h(1) = subplot(nrows,1,i);
maxn=max(n(Pl.xlims(1):Pl.xlims(2)));
Plot_n(Pl,n);
% title(['a=',num2str(Phat.a),', b=',num2str(Phat.b), ', sig=',num2str(Phat.sig), ' lam=',num2str(Phat.lam)])

% plot inferred spike trains
for r=1:1
    i=i+1; h(3+r) = subplot(nrows,1,i);
    Pl.label = 'cock'; %I{r}.P.name;
    Plot_n_MAP(Pl,I{r}.n);
    title(['a=',num2str(I{r}.P.a),', b=',num2str(I{r}.P.b), 'sig=',num2str(I{r}.P.sig),...
        ', lam=',num2str(I{r}.P.lam), ', n=',num2str(max(I{r}.n))])
end


subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)
% linkaxes(h,'x')
% linkaxes([h(end-1), h(end)])
% 
% % print fig
% wh=[7 5];   %width and height
% set(fig,'PaperPosition',[0 11-wh(2) wh]);
% print('-depsc','SimThresh')