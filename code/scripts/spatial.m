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
Nc      = 1;                        % # of cells in the ROI
neur_w  = 15;                       % width per neuron
width   = 15;                       % width of frame (pixels)
height  = Nc*neur_w;                % height of frame (pixels)
Npixs   = width*height;             % # pixels in ROI
x       = linspace(-5,5,height);
y       = linspace(-5,5,width);
[X,Y]   = meshgrid(x,y);            
g1      = zeros(Npixs,Nc);
g2      = 0*g1;
Sigma1  = diag([1,1])*2;            % var of positive gaussian
Sigma2  = diag([1,1])*2.5;          % var of negative gaussian
mu      = [0 0];                    % means of gaussians for each cell (distributed across pixel space)
w       = 1;                        % weights of each filter


ix0 = [97:99 112:114 127:129];      % typical spatial filter
for i=1:Nc                          % spatial filter
    g1(:,i)  = w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma1);
    g2(:,i)  = w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma2);
end

% 2) set simulation metadata
V.T       = 400;                    % # of time steps
V.dt      = 0.005;                  % time step size
V.MaxIter = 0;                      % # iterations of EM to estimate params
V.Np      = Npixs;                  % # of pixels in each image
V.w       = width;                  % width of frame (pixels)
V.h       = height;                 % height of frame (pixels)
V.Nc      = Nc;                     % # cells
V.plot    = 0;                      % whether to plot filter with each iteration
V.save    = 1;                                

% 3) initialize params
for i=1:V.Nc
    P.a(:,i)=g1(:,i)-1.1*g2(:,i);
end
% figure(2), imagesc(reshape(P.a,V.w,V.h)), colormap(gray), colorbar

% P.a(56,1)=P.a(56,1)-P.a(56,1)/5;
P.b     = 0;%*P.a(:,1);             % baseline is zero

P.sig   = 0.005;                    % stan dev of noise (indep for each pixel)
C_0     = 0;                        % initial calcium
tau     = round(100*rand(V.Nc,1))/100+0.05;   % decay time constant for each cell
P.gam   = 1-V.dt./tau(1:V.Nc);
P.lam   = 5;                        % rate-ish, ie, lam*dt=# spikes per second

% 3) simulate data
n=zeros(V.T,V.Nc);
C=n;
for i=1:V.Nc
    n(1,i)      = C_0;
    n(2:end,i)  = poissrnd(P.lam(i)*V.dt*ones(V.T-1,1));    % simulate spike train
    n(n>1)      = 1;
    C(:,i)      = filter(1,[1 -P.gam(i)],n(:,i));           % calcium concentration
end
Z = 0*n(:,1);
F = C*P.a' + P.b + P.sig*randn(V.T,Npixs);                  % fluorescence

%% 4) other stuff

MakMov  = 1;
% make movie of raw data
if MakMov==1
    for i=1:V.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(F(i,:),width,height),'spatial_Mov2.tif','tif','Compression','none','WriteMode',mod)
    end
end

GetROI  = 0;
fnum    = 0;

if GetROI
    figure(100); clf,imagesc(reshape(sum(g1-g2,2),width,height))
    for i=1:Nc
        [x y]   = ginput(4);
        ROWS    = [round(mean(y(1:2))) round(mean(y(3:4)))];                              % define ROI
        COLS    = [round(mean(x([1 4]))) round(mean(x(2:3)))];
        COLS1{i}=COLS;
        ROWS1{i}=ROWS;
        save('ROIs','ROWS1','COLS1')
    end
else
    load('../../data/ROIs.mat')
end


%% end-1) infer spike train using various approaches
qs=1:2;%:6;%[1 2 3];
MaxIter=10;
for q=qs
    GG=F; Tim=Sim;
    if q==1,
        Phat{q}=P;
        I{q}.label='True Filter';
    elseif q==2,
        Phat{q}=P;
        Phat{q}.a=0*P.a;
        [asort ix]=sort(P.a);
        Phat{q}.a(ix0) = 1;
        I{q}.label='Segmented Filter';
    end
    display(I{q}.label)
    [I{q}.n I{q}.P] = fast_oopsi(GG',Tim,Phat{q});
end
if V.save==1, save('../../data/spatial2'); end
%% end) plot results
% load('../../data/spatial2')
clear Pl
nrows   = 3+Nc;                                 % set number of rows
ncols   = 2;
h       = zeros(nrows,1);
Pl.xlims= [5 V.T-5];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.fs   = 13;
Pl.vs   = 2;
Pl.colors(1,:) = [0 0 0];
Pl.colors(2,:) = Pl.gray;
Pl.colors(3,:) = [.5 0 0];
Pl.Nc   = V.Nc;
Pl.XTicks=100:100:400;

% % movie slices
fnum = figure(1); clf,

movieslices=[];
Nframes=length(Pl.XTicks);
for i=1:Nframes
    temp=reshape(GG(Pl.XTicks(i),:),V.w,V.h);
    movieslices=[movieslices; temp];
end
movieslices=60*z1(movieslices);

for i=1:Nframes
    subplot('Position',[0.132+(i-1)*0.21 0.6 .2 .28])
    if i==3
        imagesc(reshape(P.a,V.w,V.h))
        colormap(gray)
        title([{'True'}; {'Spatial Filter'}],'FontSize',Pl.fs)
    elseif i==2 
        imagesc(reshape(Phat{2}.a,V.w,V.h))
        title([{'Typical'}; {'Spatial Filter'}],'FontSize',Pl.fs)
        set(gca,'YTick',[],'XTick',[])
    elseif i==4 
        imagesc(reshape(mean(F),V.w,V.h))
        title('Mean Frame','FontSize',Pl.fs)
        set(gca,'YTick',[],'XTick',[])
    elseif i==1
        imagesc(reshape(F(100,:),V.w,V.h))
        title([{'Example Frame'}],'FontSize',Pl.fs)
        set(gca,'YTick',[],'XTick',[])
    end
end

for q=qs
    i=q+3;
    if q==1, p=2; else p=1; end
    % plot fluorescence data
    i=i+1; h(i) = subplot(nrows,ncols,i);
    if q==1,
        title([{'True Projection'}],'FontSize',Pl.fs+2);
        Pl.label = 'fluorescence';
    else
        Pl.label=[];
        Pl.interp = 'none';
        title([{'Typical Projection'}],'FontSize',Pl.fs+2)
    end
    Pl.color = 'k';
    Plot_nX(Pl,(Phat{q}.a\F')');

    % plot inferred spike trains
    if q==1, 
        Pl.label = [{'fast'}; {'filter'}]; 
        Pl.interp = 'none';
    else
        Pl.label=[]; 
    end
    i=i+2; h(i) = subplot(nrows,ncols,i);
    Pl.col(2,:)=[0 0 0];
    Pl.gray=[.5 .5 .5];
    hold on
    Plot_n_MAP(Pl,I{q}.n);

    % set xlabel stuff
    subplot(nrows,ncols,i)
    set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*V.dt,'FontSize',Pl.fs)
    xlabel('Time (sec)','FontSize',Pl.fs)
end

if V.save==1 % print fig
    wh=[7 5];   %width and height
    DirName = '../../figs/';
    FileName = 'spatial2';
    PrintFig(wh,DirName,FileName);
end