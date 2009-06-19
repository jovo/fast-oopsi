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

clear, clc

% 1) generate spatial filters

% % stuff required for each spatial filter
Nc      = 1;                                % # of cells in the ROI
neur_w  = 15;                               % width per neuron
width   = 15;                               % width of frame (pixels)
height  = Nc*neur_w;                        % height of frame (pixels)
Npixs   = width*height;                     % # pixels in ROI
x       = linspace(-5,5,height);
y       = linspace(-5,5,width);
[X,Y]   = meshgrid(x,y);
g1      = zeros(Npixs,Nc);
g2      = 0*g1;
Sigma1  = diag([1,1])*2;                    % var of positive gaussian
Sigma2  = diag([1,1])*2.5;                    % var of negative gaussian
mu      = [0 0]; %[1 1]'*linspace(-2,2,Nc); % means of gaussians for each cell (distributed across pixel space)
w       = Nc:-1:1;                          % weights of each filter

% spatial filter
for i=1:Nc
    g1(:,i)  = w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma1);
    g2(:,i)  = w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma2);
end
a_b = sum(g1-g2,2);

% 2) set simulation metadata
Sim.T       = 1000;                              % # of time steps
Sim.dt      = 0.005;                            % time step size
Sim.MaxIter = 0;                                % # iterations of EM to estimate params
Sim.Np      = Npixs;                            % # of pixels in each image
Sim.w       = width;                            % width of frame (pixels)
Sim.h       = height;                           % height of frame (pixels)
Sim.Nc      = Nc;                               % # cells
Sim.plot    = 0;                                % whether to plot filter with each iteration
Sim.thresh  = 1;

% 3) initialize params
for i=1:Sim.Nc
    P.a(:,i)=g1(:,i)-1.1*g2(:,i);
end
P.b     = 0; %*P.a(:,1);                           % baseline is zero

P.sig   = 0.005;                                 % stan dev of noise (indep for each pixel)
C_0     = 0;                                    % initial calcium
tau     = .5; %round(100*rand(Sim.Nc,1))/100+0.05;   % decay time constant for each cell
P.gam   = 1-Sim.dt./tau(1:Sim.Nc);
P.lam   = 10;%round(10*rand(Sim.Nc,1))+5;           % rate-ish, ie, lam*dt=# spikes per second

% 3) simulate data
n=zeros(Sim.T,Sim.Nc);
C=n;
for i=1:Sim.Nc
    n(1,i)      = C_0;
    n(2:end,i)  = poissrnd(P.lam(i)*Sim.dt*ones(Sim.T-1,1));    % simulate spike train
    n(n>1)      = 1;
    C(:,i)      = filter(1,[1 -P.gam(i)],n(:,i));               % calcium concentration
end
Z = 0*n(:,1);
F = C*P.a' + P.b + P.sig*randn(Sim.T,Npixs);               % fluorescence

 imagesc(reshape(F(500,:),Sim.w,Sim.h))
 
%% 4) other stuff

MakMov  = 1;
% make movie of raw data
if MakMov==1
    for i=1:Sim.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(F(i,:),width,height),'spatial_Mov.tif','tif','Compression','none','WriteMode',mod)
    end
end

%% 5) infer spike train using various approaches
qs=1:2;%:6;%[1 2 3];
MaxIter=10;
for q=qs
    GG=F; Tim=Sim; Phat{q}=P;
    if q==1,
        I{q}.label='True Parameters';
        [I{q}.n I{q}.P] = FOOPSI_v3_05_01(GG',Phat{q},Tim);
    elseif q==2,
        I{q}.label='Estimated Parameters';
        Est.Thresh  = 1;
        Est.a       = 1;
        Est.b       = 1;
        Est.sig     = 1;
        Est.lam     = 1;
        Tim.MaxIter = 10;
        Est.Plot    = 1;
        PP          = P;
        PP.b        = 1;
%         PP.gam      = 2*P.gam;
        PP.lam      = 2*P.lam;
        PP.sig      = 2*P.sig;
        [I{q}.n I{q}.P] = FOOPSI_v3_05_01(GG',PP,Tim,Est);
    end
    display(I{q}.label)
    %     [I{q}.n I{q}.P] = FOOPSI2_59(GG,Phat{q},Tim);
end

%% end) plot results
clear Pl
nrows   = 3;                                  % set number of rows
ncols   = numel(qs);
h       = zeros(nrows,1);
Pl.xlims= [1 201]+100;                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 2;
Pl.colors(1,:) = [0 0 0];
Pl.colors(2,:) = Pl.gray;
Pl.colors(3,:) = [.5 0 0];
Pl.Nc   = Sim.Nc;
Pl.XTicks=[200 400 600];

figure(3), clf
for q=qs

    % plot spatial filter
    i=q; h(i) = subplot(nrows,ncols,i);
    imagesc(reshape(Phat{q}.a,Sim.w,Sim.h)),
    title(I{q}.label)
    if q==1,
        ylab=ylabel([{'Spatial'}; {'Filter'}]);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    else
        set(gca,'XTick',[],'YTick',[])
    end
    colormap('gray')

    % plot fluorescence data
    i=i+ncols; h(i) = subplot(nrows,ncols,i);
    Pl.color = 'k';
    if q==1, 
        Pl.label='$F^{proj}_{1:T}$';        
    else
        Pl.label=[]; %'$F^{est}_{1:T}$'; 
    end
    plot(z1(Phat{q}.a\F(Pl.xlims(1):Pl.xlims(2),:)'),'k','LineWidth',Pl.lw)
    axis([0 Pl.xlims(2)-Pl.xlims(1) 0 1])
    box off
    set(gca,'XTick',[0 100 200],'XTickLabel',[]); %[0 100 200]*Sim.dt)
    set(gca,'YTick',[0 .5 1],'YTickLabel',[])
    
    ylab=ylabel(Pl.label,'Interpreter','latex','FontSize',Pl.fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

    %     Plot_nX(Pl,(Phat{q}.a\F')');

    % plot inferred spike trains
    if q==1, Pl.label = [{'$\mathbf{n}^{FAND}$'}];
    else Pl.label=[]; end
    Pl.interp='latex';
    i=i+ncols; h(i) = subplot(nrows,ncols,i);
    Pl.col(2,:)=[0 0 0];
    Pl.gray=[.5 .5 .5];
    hold on
    Plot_n_MAP(Pl,I{q}.n);

    % set xlabel stuff
    subplot(nrows,ncols,i)
    set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
    xlabel('Time (sec)','FontSize',Pl.fs)
    %     linkaxes(h,'x')

end

% print fig
wh=[7 5];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
FigName = '../../figs/spatial_EM';
print('-depsc',FigName)
print('-dpdf',FigName)