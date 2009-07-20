clear, clc, %close all
datapath = '/Users/joshyv/Research/oopsi/pop-oopsi/data/';
dataset = 2;
switch dataset
    case 1
        Im.path = [datapath 'rafa/tanya/051809/s1m1/'];
        Im.fname= 's1m1';
        Im.Fura = 1;
    case 2
        Im.fname= '20081126_13_05_43_orientation_Bruno_reg1_ori1_135umdepth';
        Im.matname = [datapath Im.fname];
    case 3
        Im.path = [datapath 'mrsic-flogel/'];
        Im.fname= '20081126_13_15_31_natural_reg1_nat_135umdepth';
        Im.Fura = 0;
end
load(Im.matname)

%% 2) set simulation metadata
Meta.MaxIter= 25;                               % # iterations of EM to estimate params
Meta.Np     = 1;                     % # of pixels in each image
% Meta.matname= Im.matname;
Meta.dt     = 1/15;

User.Nc     = 1;                               % # cells
User.Plot   = 1;                               % whether to plot filter with each iteration
User.Thresh = 1;
User.lam    = 1;
           
for k=1:numel(F)

    %% set parameters
    
    I{k}.P0.a   = mean(F{k});
    I{k}.P0.b   = I{k}.P0.a/10; 
    I{k}.P0.sig = std(F{k});
    I{k}.P0.lam = 1/30;
    tau         = 1;
    I{k}.P0.gam = (1-Meta.dt/tau);

    %% infer spike train using various approaches
    [I{k}.n I{k}.P] = FOOPSI_v3_05_01(F{k},I{k}.P0,Meta,User);

end

%% 5) plot results
clear Pl
ncols   = 3;                                  % set number of rows
nrows   = 4; %numel(qs);
h       = zeros(nrows,1);
Pl.xlims= [5 Meta.T-30];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
n       = zeros(Meta.T,1); n(Meta.spt)=1;
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 4;
Pl.colors(1,:) = [0 0 0];
Pl.colors(2,:) = Pl.gray;
Pl.colors(3,:) = [.5 0 0];
Pl.Nc   = User.Nc;
Pl.XTicks=[100:100:Meta.T];
Pl.interp = 'latex';

fnum=figure(1); clf

MF=60*z1(MeanFrame);
PA=(I{q}.P.a)/max(I{q}.P.a);

subplot('Position',[.05 .44 .9 .5]);
imagesc(reshape(MF,MovInf(1).Height,MovInf(1).Width))
title('Mean Movie Frame','FontSize',Pl.fs)

for q=qs

    % plot spatial filter
    i=7; %h(i) =
    subplot('Position',[0.05 0.05 .35 .3]) %nrows,ncols,[7]);
    imagesc(reshape(PA,Meta.h,Meta.w)'),
    colormap('gray')
    title('Spatial Filter','FontSize',Pl.fs)

    % plot fluorescence data
    subplot('Position',[.45 .215 .5 .12])
    Pl.color = 'k';
    %         Pl.label=[{'Fluorescence'}; {'Projection'}];
    Pl.label=['$\widetilde{F}$'];
    Plot_nX(Pl,F*Phat{q}.a);
    set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)

    % plot inferred spike trains
    subplot('Position',[.45 .09 .5 .12])
    Pl.col(2,:)=[0 0 0];
    Pl.gray=[.5 .5 .5];
    hold on
    Pl.label=['$\widehat{n}$'];
    Plot_n_MAP(Pl,I{q}.n);
    set(gca,'XTick',Pl.XTicks,'XTickLabel',(Pl.XTicks)*Meta.dt,'FontSize',Pl.fs)
    xlabel('Time (sec)','FontSize',Pl.fs)
end

% print fig
wh=[7 7];   %width and height
DirName = '../graphics/';
FileName = 'spatial_data';
PrintFig(wh,DirName,FileName);
