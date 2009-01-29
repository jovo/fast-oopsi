% this script generates a simulation of a movie containing a single cell
% using the following generative model:
% 
% F_t = a*C_t + b + sig*eps_t, eps_t ~ N(0,I)     
% C_t = gam*C_{t-1} + n_t,      n_t ~ Poisson(lam*dt)
% 
% where a,b,I are p-by-q matrices.  
% we let b=0 and a be the difference of gaussians (yielding a zero mean
% matrix)

clear, clc

% 1) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.MaxIter = 5;

% 2) initialize spatial filter
Sim.Nrows = 30;
Sim.Ncols = Sim.Nrows;
x1      = linspace(-5,5,Sim.Nrows);
x2      = x1;
[X1,X2] = meshgrid(x1,x2);

mu      = [0 0];
Sigma   = diag([1,1])*1;
a       = mvnpdf([X1(:) X2(:)],mu,Sigma);

mu      = [0 0];
Sigma   = diag([1,1])*2;
b       = mvnpdf([X1(:) X2(:)],mu,Sigma);

figure(1), clf, imagesc(reshape(a-b,Sim.Nrows,Sim.Ncols));

% 3) initialize other params
P.a     = a-b;
P.b     = 0*P.a;
P.sig   = .005;                                  % stan dev of noise
C_0     = 0;
tau     = 0.05;                                 % decay time constant
P.gam   = 1-Sim.dt/tau;
P.lam   = 10;                                   % rate-ish, ie, lam*dt=# spikes per second

% 3) simulate data
n = poissrnd(P.lam*Sim.dt*ones(Sim.T-1,1));     % simulate spike train
% n(n>1)=1;
n = [C_0; n];                                   % set initial calcium
C = filter(1,[1 -P.gam],n);                     % calcium concentration
Z = 0*C;
F = C*P.a'+Z*P.b'+P.sig*randn(Sim.T,length(P.a));             % fluorescence

% 4) estimate spatial filter from real spikes
X       = [C 1+Z];
Phat    = P;
for ii=1:numel(P.a)
    Y   = F(:,ii);
    B   = X\Y;
    Phat.a(ii) = B(1);
    Phat.b(ii) = B(2);
end

% for i=1:Sim.T
%     imwrite(reshape(F(i,:),Sim.Nrows,Sim.Nrows)',['Sim1Cell.tif'],'tif','Compression','none','WriteMode','append')
% end

% 4) infer spike train
for q=1:5
    if q==1
        F2=mean(F,2);
        P2=P;
        P2.a=mean(P.a);
        P2.b=mean(P.b);
        [I{q}.n I{q}.P] = FOOPSI2_57(F2,P2,Sim);
        label = [{'Uniform'}; {'Filter'}];
    elseif q==2
        [I{q}.n I{q}.P] = FOOPSI2_57(F,P,Sim);
        label = [{'True'}; {'Filter'}];
    elseif q==3
        [I{q}.n I{q}.P] = FOOPSI2_57(F,Phat,Sim);
        label = [{'Spike'}; {'Filter'}];
    elseif q==4
        MeanFrame=mean(F)';
        P3=P;
        P3.a=MeanFrame; %reshape(MeanFrame,Sim.Nrows,Sim.Ncols);
        [I{q}.n I{q}.P] = FOOPSI2_57(F,P3,Sim);
        label = [{'Mean'}; {'Filter'}];
    else    
        Sim.MaxIter=5;
        [I{q}.n I{q}.P] = FOOPSI2_57(F,P3,Sim);
        label = [{'Fluor'}; {'Filter'}];        
    end
        I{q}.P.label = label;
end

%%
% 5) plot results
fig     = figure(2); clf,
nrows   = 3+q;                                  % set number of rows
h       = zeros(nrows,1);
Pl.xlims= [5 Sim.T];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 2;

% plot fluorescence data
i=1; h(1) = subplot(nrows,1,i);
Pl.label = 'Fluorescence';
Pl.color = 'k';
Plot_X(Pl,F*mean(F)');

% plot calcium
i=i+1; h(2) = subplot(nrows,1,i);
Pl.label = 'Calcium';
Pl.color = Pl.gray;
Plot_X(Pl,C);

% plot spike train
i=i+1; h(3) = subplot(nrows,1,i);
maxn=max(n(Pl.xlims(1):Pl.xlims(2)));
Plot_n(Pl,n);
% title(['a=',num2str(Phat.a),', b=',num2str(Phat.b), ', sig=',num2str(Phat.sig), ' lam=',num2str(Phat.lam)])

% plot inferred spike trains
for r=1:q
    i=i+1; h(3+r) = subplot(nrows,1,i);
    Pl.label = I{r}.P.label;
    Plot_n_MAP(Pl,I{r}.n);
    %     title(['a=',num2str(mean(I{r}.P.a)),', b=',num2str(I{r}.P.b), ', sig=',num2str(I{r}.P.sig),...
    %         ', lam=',num2str(I{r}.P.lam), ', n_{max}=',num2str(max(I{r}.n))])
end


subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)
linkaxes(h,'x')

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','Spatial')

%%
fig = figure(3); clf
nrows = 3;
ncols = 3;

mx = max([MeanFrame; P.a; P.b; Phat.a; Phat.b; I{q}.P.a; I{q}.P.b]);
mn = min([MeanFrame; P.a; P.b; Phat.a; Phat.b; I{q}.P.a; I{q}.P.b]);

subplot(nrows,ncols,1), imagesc(reshape((MeanFrame-mn)/mx,Sim.Nrows,Sim.Ncols))%, colorbar
subplot(nrows,ncols,2), imagesc(reshape((P.a-mn)/mx,Sim.Nrows,Sim.Ncols))%, colorbar
subplot(nrows,ncols,3), imagesc(reshape((P.b-mn)/mx,Sim.Nrows,Sim.Ncols))

subplot(nrows,ncols,4), imagesc(reshape((MeanFrame-mn)/mx,Sim.Nrows,Sim.Ncols))%, colorbar
subplot(nrows,ncols,5), imagesc(reshape((Phat.a-mn)/mx,Sim.Nrows,Sim.Ncols))%, colorbar
subplot(nrows,ncols,6), imagesc(reshape((Phat.b-mn)/mx,Sim.Nrows,Sim.Ncols))

subplot(nrows,ncols,7), imagesc(reshape((MeanFrame-mn)/mx,Sim.Nrows,Sim.Ncols))%, colorbar
subplot(nrows,ncols,8), imagesc(reshape((I{q}.P.a-mn)/mx,Sim.Nrows,Sim.Ncols))%, colorbar
subplot(nrows,ncols,9), imagesc(reshape((I{q}.P.b-mn)/mx,Sim.Nrows,Sim.Ncols))
