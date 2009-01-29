% this script generates a simulation of a movie containing a single cell
% using the following generative model:
%
% F_t = \sum_i a_i*C_{i,t} + b + sig*eps_t, eps_t ~ N(0,I)
% C_{i,t} = gam*C_{i,t-1} + n_{i,t},      n_{i,t} ~ Poisson(lam_i*dt)
%
% where ai,b,I are p-by-q matrices.
% we let b=0 and ai be the difference of gaussians (yielding a zero mean
% matrix)
%
% Note that this version is can generate up to 2 cells (because we manually
% assign each cell a distribution

clear, clc

MakMov = 0;
% 1) generate spatial filters

% stuff required for each spatial filter
Nc      = 3;                                % # of cells
Nrows   = 20;                               % # rows in each image frame
colpern = 13;                               % # columns per neuron
Ncols   = Nc*colpern;                       % # cols in each image frame
Npixs   = Nrows*Ncols;
x1      = linspace(-5,5,Ncols);
x2      = linspace(-5,5,Nrows);
[X1,X2] = meshgrid(x1,x2);
a       = zeros(Npixs,Nc);
b       = 0*a;
Sigma1  = diag([1,1])*1;
Sigma2  = diag([1,1])*2;
mu     = [1 1]'*linspace(-4,4,Nc);

% spatial filter
for i=1:Nc
    %     mu      = randn(1,2)*2 %
    %     mu      = rand(1,2).*[i*colpern Nrows] + [(i-1)*colpern 0]
    a(:,i)  = mvnpdf([X1(:) X2(:)],mu(:,i)',Sigma1);
    b(:,i)  = mvnpdf([X1(:) X2(:)],mu(:,i)',Sigma2);
end

% 2) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.MaxIter = 0;
Sim.Np      = Npixs;                            % # of pixels in each image
Sim.Nc      = Nc;                               % # cells

% 3) initialize params
P.a     = zeros(Sim.Np,Sim.Nc);
for i=1:Sim.Nc
    P.a(:,i)=a(:,i)-b(:,i);
end
P.b     = 0*P.a(:,1);

P.sig   = 0.05;                                  % stan dev of noise
C_0     = 0;
tau     = round(100*rand(Sim.Nc,1))/100+0.05;                                 % decay time constant
P.gam   = 1-Sim.dt./tau(1:Sim.Nc);
P.lam   = round(10*rand(Sim.Nc,1))+5; %10*ones(Sim.Nc,1);                                   % rate-ish, ie, lam*dt=# spikes per second

% 3) simulate data
n=zeros(Sim.T,Sim.Nc);
C=n;
for i=1:Sim.Nc
    n(1,i)      = C_0;
    n(2:end,i)  = poissrnd(P.lam(i)*Sim.dt*ones(Sim.T-1,1));    % simulate spike train
    C(:,i)      = filter(1,[1 -P.gam(i)],n(:,i));               % calcium concentration
end
Z = 0*n(:,1);
F = C*P.a' + (1+Z)*P.b'+P.sig*randn(Sim.T,Npixs);                   % F = C1*P.a1'+ C2*P.a2' + Z*P.b'+P.sig*randn(Sim.T,Npixs);             % fluorescence

figure(2), clf, bar(n(:,1),'k'),
if Nc>1, hold on, bar(-n(:,2),'b'), end

% 4) estimate spatial filter from real spikes
X       = [C 1+Z];
Phat{1} = P;
for i=1:Npixs
    Y   = F(:,i);
    B   = X\Y;
    for j=1:Nc
        Phat{1}.a(i,j) = B(j);
    end
    Phat{1}.b(i) = B(end);
end
fig=figure(1); clf, ncols=2+Nc;
mn = min([min(sum(a-b,2)) min(Phat{1}.a(:)) min(Phat{1}.b(:))]);
mx = (max([max(sum(a-b,2)) max(Phat{1}.a(:)) max(Phat{1}.b(:))])-mn)/60;

subplot(1,ncols,1), image(reshape((sum(a-b,2)-mn)/mx,Nrows,Ncols))%, colorbar
for i=1:Nc
    subplot(1,ncols,1+i), image(reshape((Phat{1}.a(:,i)-mn)/mx,Nrows,Ncols))%, colorbar
    title(['a' num2str(i)])
end
subplot(1,ncols,ncols), image(reshape((Phat{1}.b-mn)/mx,Nrows,Ncols)), title('b')%, colorbar
wh=[7 3];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-deps','Multi_Filters')

%% 5) try denoising using SVD of whole movie

% make movie of raw data
if MakMov==1
    for i=1:Sim.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(F(i,:),Nrows,Ncols),'Multi_Mov.tif','tif','Compression','none','WriteMode',mod)
    end
end

[U,S,V] = svd(F,0);
PCs = 1:3;
Denoised = U(:,PCs)*S(PCs,PCs)*V(:,PCs)';

if MakMov==1
    for i=1:Sim.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(Denoised(i,:),Nrows,Ncols),'Multi_Mov_Denoised.tif','tif','Compression','none','WriteMode',mod)
    end
end

mn = min([min(sum(a-b,2)) min(min(V(:,1:Nc)))]);
mx = (max([max(sum(a-b,2)) max(max(V(:,1:Nc)))])-mn)/60;

figure(3), clf,ncols=2+Nc;
subplot(1,ncols,1), image(reshape((sum(a-b,2)-mn)/mx,Nrows,Ncols))%, colorbar
Phat{2}=P;
for i=1:Nc
    Phat{2}.a(:,i) = V(:,i);
    subplot(1,ncols,1+i), image(reshape((Phat{2}.a(:,i)-mn)/mx,Nrows,Ncols))%, colorbar
end

%% 6) try denoising by doing SVD on each ROI

Phat{3}=P;
r=1:round(Nrows/Nc):Nrows;
c=1:round(Ncols/Nc):Ncols;
for i=1%:Nc
    % get indices
    ROWS    = [r(i) r(i)+5];
    COLS    = [c(i) c(i)+13];
    NROWS   = diff(ROWS)+1;
    NCOLS   = diff(COLS)+1;
    ind=1:NROWS;
    for j=1:NCOLS-1
        ind=[ind j*Nrows+1:j*Nrows+NROWS];
    end

    FF          = F(:,ind);
    [U,S,V]     = svd(FF,0);
    Phat{3}.a(:,i)   = 0;
    Phat{3}.a(ind,i) = V(:,1);
end

mn = min([min(sum(a-b,2)) min(P.a(:)) min(Phat{3}.a(:))]);
mx = (max([max(sum(a-b,2)) max(P.a(:)) max(Phat{3}.a(:))])-mn)/60;

figure(4), clf,ncols=1+Nc;
subplot(1,ncols,1), image(reshape((sum(a-b,2)-mn)/mx,Nrows,Ncols))%, colorbar
for i=1:Nc
    subplot(1,ncols,1+i), image(reshape((Phat{3}.a(:,i)-mn)/mx,Nrows,Ncols))%, colorbar
end

%% 4) infer spike train

for q=1:3
    if q==1
        [I{q}.n I{q}.P] = FOOPSI2_59(F,Phat{1},Sim);
        I{q}.P.label    = 'True Filter';
    elseif q==2
        [I{q}.n I{q}.P] = FOOPSI2_59(F,Phat{2},Sim);
        I{q}.P.label    = 'SVD Filter';
    elseif q==3
        [I{q}.n I{q}.P] = FOOPSI2_59(F,Phat{3},Sim);
        I{q}.P.label    = 'ROI SVD Filter';
    end
end

%%
% 5) plot results
clear Pl
nrows   = 2+Nc;                                  % set number of rows
h       = zeros(nrows,1);
Pl.xlims= [5 Sim.T];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 2;
Pl.colors(1,:) = [0 0 0];
Pl.colors(2,:) = Pl.gray;
Pl.colors(3,:) = [.5 0 0];
Pl.Nc   = Sim.Nc;

for r=1:q
    fig = figure(4+r); clf,

    % plot fluorescence data
    i=1; h(i) = subplot(nrows,1,i);
    Pl.label = [{'Spatially'}; {'Filtered'}; {'Fluorescence'}];
    Pl.color = 'k';
    Plot_nX(Pl,F*Phat{r}.a);
    title(I{r}.P.label)

    % plot calcium
    i=i+1; h(i) = subplot(nrows,1,i);
    Pl.label = 'Calcium';
    Pl.color = Pl.gray;
    Plot_nX(Pl,C);

    % plot inferred spike trains
    Pl.label = [{'Spike'}; {'Train'}];
    for j=1:Nc
        i=i+1; h(i) = subplot(nrows,1,i);
        Pl.j=j;
        Plot_n_MAPs(Pl,I{r}.n(:,j));
    end

    % set xlabel stuff
    subplot(nrows,1,nrows)
    set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
    xlabel('Time (sec)','FontSize',Pl.fs)
    linkaxes(h,'x')

    % print fig
    wh=[7 5];   %width and height
    set(fig,'PaperPosition',[0 11-wh(2) wh]);
    print('-depsc',['Multi_Spikes' num2str(r)])
end