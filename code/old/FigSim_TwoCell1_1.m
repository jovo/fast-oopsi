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

% 1) generate spatial filters

% stuff required for each spatial filter
Nc      = 2;                                % # of cells
Nrows   = 30;
Ncols   = 30;
Npixs   = Nrows*Ncols;
x1      = linspace(-5,5,Nrows);
x2      = x1;
[X1,X2] = meshgrid(x1,x2);
a       = zeros(Npixs,Nc);

% spatial filter 1
mu      = [-2 -2];
Sigma   = diag([1,1])*1;
a(:,1)  = mvnpdf([X1(:) X2(:)],mu,Sigma);
Sigma   = diag([1,1])*2;
b(:,1)  = mvnpdf([X1(:) X2(:)],mu,Sigma);

% spatial filter 2
mu      = [2 2];
Sigma   = diag([1,1])*1;
a(:,2)  = mvnpdf([X1(:) X2(:)],mu,Sigma);
Sigma   = diag([1,1])*2;
b(:,2)  = mvnpdf([X1(:) X2(:)],mu,Sigma);

figure(1), clf, ncols=2+Nc;
subplot(1,ncols,1), imagesc(reshape(sum(a-b,2),Nrows,Ncols))%, colorbar

% 2) set simulation metadata
Sim.T       = 500;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.MaxIter = 0;
Sim.Np      = Npixs;                            % # of pixels in each image
Sim.Nc      = Nc;                               % # cells

% 3) initialize params
P.a     = zeros(Sim.Np,Sim.Nc);                  % P.a1    = a1-b1; % P.a2    = a2-b2; % [a1-b1 a2-b2];
for i=1:Sim.Nc
    P.a(:,i)=a(:,i)-b(:,i);
end
P.b     = 0*P.a(:,1);

P.sig   = 0.05;                                  % stan dev of noise
C_0     = 0;
tau     = [0.05; 0.5];                                 % decay time constant
P.gam   = 1-Sim.dt./tau(1:Sim.Nc);
P.lam   = 10*ones(Sim.Nc,1);                                   % rate-ish, ie, lam*dt=# spikes per second

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
Phat    = P;
for i=1:Npixs
    Y   = F(:,i);
    B   = X\Y;
    for j=1:Nc
        Phat.a(i,j) = B(j);
    end
    Phat.b(i) = B(end);
end
fig=figure(1);
for i=1:Nc
    subplot(1,ncols,1+i), imagesc(reshape(Phat.a(:,i),Nrows,Ncols))%, colorbar
    title(['a' num2str(i)])
end
subplot(1,ncols,ncols), imagesc(reshape(Phat.b,Nrows,Ncols)), title('b')%, colorbar
wh=[7 3];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
% print('-dpdf','svd_good')

%% 5) try denoising using SVD
%
% % make movie of raw data
% for i=1:Sim.T
%     if i==1, mod='overwrite'; else mod='append'; end
%     imwrite(reshape(F(i,:),Nrows,Nrows)',['Sim2Cell.tif'],'tif','Compression','none','WriteMode',mod)
% end
%
% [U,S,V] = svd(F,0);
% PCs = 1:2;
% Denoised = U(:,PCs)*S(PCs,PCs)*V(:,PCs)';
%
% for i=1:Sim.T
%     if i==1, mod='overwrite'; else mod='append'; end
%     imwrite(reshape(Denoised(i,:),Nrows,Nrows)',['Sim2Cell_Denoised.tif'],'tif','Compression','none','WriteMode',mod)
% end
%
% figure(3), clf,
% subplot(1,4,1), imagesc(reshape(P.a1+P.a2,Nrows,Ncols))%, colorbar
% subplot(1,4,2), imagesc(reshape(V(:,1),Nrows,Ncols))%, colorbar
% subplot(1,4,3), imagesc(reshape(V(:,2),Nrows,Ncols))%, colorbar
% subplot(1,4,4), imagesc(reshape(V(:,3),Nrows,Ncols))%, colorbar

%% 4) infer spike train
q=1;
[I{q}.n I{q}.P] = FOOPSI2_59(F,P,Sim);
I{q}.P.label    = 'ass';
% for q=1:5
%     if q==1
%         F2=mean(F,2);
%         P2=P;
%         P2.a=mean(P.a);
%         P2.b=mean(P.b);
%         [I{q}.n I{q}.P] = FOOPSI2_57(F2,P2,Sim);
%         label = [{'Uniform'}; {'Filter'}]
%     elseif q==2
%         [I{q}.n I{q}.P] = FOOPSI2_57(F,P,Sim);
%         label = [{'True'}; {'Filter'}]
%     elseif q==3
%         [I{q}.n I{q}.P] = FOOPSI2_57(F,Phat,Sim);
%         label = [{'Spike'}; {'Filter'}]
%     elseif q==4
%         MeanFrame=mean(F)';
%         P3=P;
%         P3.a=MeanFrame; %reshape(MeanFrame,Nrows,Ncols);
%         [I{q}.n I{q}.P] = FOOPSI2_57(F,P3,Sim);
%         label = [{'Mean'}; {'Filter'}]
%     else
%         Sim.MaxIter=5;
%         [I{q}.n I{q}.P] = FOOPSI2_57(F,P3,Sim);
%         label = [{'Fluor'}; {'Filter'}]
%     end
%         I{q}.P.label = label;
% end

%%
% 5) plot results
clear Pl
fig     = figure(2); clf,
nrows   = 3+Nc;                                  % set number of rows
h       = zeros(nrows,1);
Pl.xlims= [5 Sim.T];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 2;
Pl.colors(1,:) = [0 0 0];
Pl.colors(2,:) = Pl.gray;
Pl.Nc   = Sim.Nc;

% plot fluorescence data
i=1; h(i) = subplot(nrows,1,i);
Pl.label = 'Fluorescence'; %[{'Optimally'}; {'Filtered'}; {'Fluorescence'}];
Pl.color = 'k';
Plot_nX(Pl,F*P.a);

% plot calcium
i=i+1; h(i) = subplot(nrows,1,i);
Pl.label = 'Calcium';
Pl.color = Pl.gray;
Plot_nX(Pl,C);

% plot spike train
i=i+1; h(i) = subplot(nrows,1,i);
Plot_2n(Pl,n);

% plot inferred spike trains
for j=1:Nc
    i=i+1; h(i) = subplot(nrows,1,i);
    cla, hold on
    bar(n(:,j))
    plot(z1(I{q}.n(:,j)))
    axis('tight')
end

% i=i+1; h(i) = subplot(nrows,1,i);
% cla, hold on
% bar(n(:,2))
% plot(z1(I{q}.n(:,2)))
% axis('tight')

% set xlabel stuff
subplot(nrows,1,nrows)
set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
xlabel('Time (sec)','FontSize',Pl.fs)
linkaxes(h,'x')

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','TwoSpatialCells')