% this script is a wrapper that loads some data (or simulates it) and then
% calls the most recent FOOPSI function to infer spike trains and
% parameters

clear, clc,
fname = 'spatial_multi2';

%% set parameters

% generate spatial filters
Nc      = 2;                                % # of cells in the ROI
neur_w  = 8;                               % height per neuron
height  = 10;                               % height of frame (pixels)
width   = Nc*neur_w;                        % width of frame (pixels)
Npixs   = height*width;                     % # pixels in ROI
x       = linspace(-5,5,width);
y       = linspace(-5,5,height);
[X,Y]   = meshgrid(x,y);
g1      = zeros(Npixs,Nc);
g2      = 0*g1;
Sigma1  = diag([1,1])*1.5;                    % var of positive gaussian
Sigma2  = diag([1,1])*5;                    % var of negative gaussian
sp      = 2.0;
mu      = [1 1]'*linspace(-sp,sp,Nc);     % means of gaussians for each cell (distributed across pixel space)
w       = ones(1,Nc); %linspace(1,Nc,Nc); %Nc:-1:1;             % weights of each filter
for i=1:Nc
    g1(:,i)  = w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma1);
    g2(:,i)  = 0*w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma2);
end

for i=1:Nc
    P.a(:,i)=g1(:,i)-g2(:,i);
end
P.b     = 0.01*mean(P.a); % baseline is a scaled down version of the sum of the spatial filters

% set meta parameters
V.T      = 4800;                         % # of time steps
V.dt     = 1/60;                           % time step size
V.Npixels= Npixs;                        % # of pixels in each image
V.Ncells = Nc;                           % # cells per ROI
V.Poiss  = 0;                            % whether observations are poisson or gaussian
V.fast_iter_max= 25;                           % # iterations of EM to estimate params

% initialize other parameters
P.sig   = 2.5*max(P.a(:));                             % stan dev of noise (indep for each pixel)
C_0     = 0;                                % initial calcium
tau     = 0.5*ones(Nc,1);   % decay time constant for each cell
P.gam   = 1-V.dt./tau(1:Nc);             % set gam
lam     = linspace(0,4*pi,V.T);
P.lam = repmat(exp(1.3*(sin(lam))'),1,V.Ncells); % rate-ish, ie, lam*dt=# spikes per second
for j=2:V.Ncells
    P.lam(:,j) = exp(1.6*(sin(lam+pi))');
end
P.lam   = P.lam-repmat(min(P.lam),V.T,1);

%% simulate data

n=zeros(V.T,V.Ncells);                         % pre-allocate memory for spike train
C=n;                                        % pre-allocate memory for calcium
for i=1:V.Ncells
    n(1,i)      = C_0;
    n(2:end,i)  = poissrnd(P.lam(2:end,i)*V.dt);  % simulate spike train
    n(n>1)      = 1;                                            % make only 1 spike per bin
    C(:,i)      = filter(1,[1 -P.gam(i)],n(:,i));               % calcium concentration
end
F = P.a*(C+repmat(P.b,V.T,1))'+P.sig*rand(Npixs,V.T);

figure(2), clf, hold off, plot(P.a(:,1)\F), hold all, plot(P.a(:,2)\F);

save(['../../data/' fname '.mat'],'F','n','P','V')

%% generate tif

MakMov  = 0;
if MakMov==1
    FF=uint8(floor(255*z1(F)));
    for t=1:500
        if t==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(FF(:,t),height,width),['../../data/' fname '.tif'],'tif','Compression','none','WriteMode',mod)
    end
end

%% infer spike trains and parameters

% initialize parameters for estimating parameters (these are only used if MaxIter>1)
V.est_sig     = 0;
V.est_lam     = 0;
V.est_gam     = 0;
V.est_b       = 1;
V.est_a       = 1;

% initialize parameters for plotting results after each pseudo-EM step
V.fast_thresh  = 1;                            % whether to threshold spike train before estimating 'a' and 'b' (we always keep this on)
V.fast_plot    = 1;                            % whether to plot filter with each iteration
V.n       = n;                            % keep true spike times to also plot for comparison purposes
V.h       = height;                       % height of frame (pixels)
V.w       = width;                        % width of frame (pixels)


% infer spike trains using a variety of techniques
q=0;
exps=[1 2.5]; %[1 2.5];% 3 3.5];
for qq=exps
    disp(qq)
    if qq==1;                                % use true params
        PP=P;
        PP.lam=mean(PP.lam);
        FF=F;
        V.fast_iter_max = 1;
    elseif qq==2                             % use svd spatial filter
        V.fast_iter_max = 1;
        PP=P;
        PP.lam=mean(PP.lam);
        FF=F;
        [U,S,VV]=pca_approx(F',V.Ncells);
        for j=1:V.Ncells, PP.a(:,j)=VV(:,j); end
    elseif qq==2.5
        V.fast_iter_max = 10;
        PP=P;
        PP.lam=mean(PP.lam);
        FF=F;
        [U,S,VV]=pca_approx(F',V.Ncells);
        for j=1:V.Ncells, PP.a(:,j)=VV(:,j); end
    elseif qq==3                           % estimate params
        V.fast_iter_max = 10;
        PP=P;
        FF=F;%-repmat(mean(F,2),1,V.T);
        [U,S,VV]=pca_approx(FF',V.Ncells);
        for j=1:V.Ncells, PP.a(:,j)=VV(:,j); end
    elseif qq==3.5                           % estimate params
        V.fast_iter_max = 10;
        PP=P;
        FF=F-repmat(mean(F,2),1,V.T);
        [U,S,VV]=pca_approx(FF',V.Ncells);
        for j=1:V.Ncells, PP.a(:,j)=VV(:,j); end
        V.est_b=0;
    elseif qq==4                           % estimate params
        V.fast_iter_max = 25;
        PP=P;
        FF=F;
        [U,S,VV]=pca_approx(FF',V.Ncells);
        for j=1:V.Ncells, PP.a(:,j)=VV(:,j); end
        V.Plot = 1;
    elseif qq==5                           % estimate params
        V.fast_iter_max = 20;
        V.fast_thresh  = 0;
        PP=P;
        [U,S,V]=pca_approx(F',V.Ncells);
        for j=1:V.Ncells, PP.a(:,j)=V(:,j); end
    end
    q=q+1;
    [I{q}.n I{q}.P] = fast_oopsi(FF,V,PP);
end

save(['../../data/' fname],'I','V')
sound(10*sin(linspace(0,180*pi,2000)))

%% plot results
fname='spatial_multi2';
load(['../../data/' fname])
% I{3}=I{2};

Pl.g    = 0.5*ones(1,3);       % gray color
Pl.fs   = 8;                   % font size
Pl.w1   = 0.28;                 % width of subplot
Pl.wh   = [Pl.w1 Pl.w1];        % width and height of subplots
Pl.b1   = 0.67;                 % bottom of topmost subplt
Pl.l1   = 0.41;                 % left side of left subplot on right half
Pl.l2   = Pl.l1+0.3;            % left side of right subplot on right half
Pl.s1   = Pl.w1+0.043;          % space between subplots on right side
Pl.s2   = .33;                  % space between subplots on left side
Pl.ms   = 2;                    % marker size
Pl.lw   = 2;                    % line width
Pl.n    = n; Pl.n(Pl.n==0)=NaN; % true spike train (0's are NaN's so they don't plot)
Pl.T    = V.T;
Pl.c    = [0 0 0; Pl.g; 1 1 1]; % colors: black, grey, white
Pl.m    = ['v','v'];
Pl.xlim = [2000 2500]-1000;        % limits of x-axis
Pl.shift= [0 .07];
Pl.x_range = Pl.xlim(1):Pl.xlim(2);
Pl.XTick = [Pl.xlim(1):2/V.dt:Pl.xlim(2)];
Pl.XTickLabel = Pl.XTick;

% show how our estimation procedure given truth and when estimating spatial filter

figure(1), clf, hold on
nrows   = 2; 
ncols   = 5; 

Fmean=mean(F);

% make images on same scale
immax=max((P.a(:)));
immin=min((P.a(:)));
for q=1:2
    immax=max(immax,max(I{q}.P.a(:)));
    immin=min(immin,min(I{q}.P.a(:)));
end

PP.a=60*(P.a-immin)/(immax-immin);
for q=1:2
    for j=1:2
        EI{q}.P.a(:,j)=60*(I{q}.P.a(:,j)-immin)/(immax-immin);
    end
end

% plot mean image frame
subplot(nrows,ncols,1)
imagesc(reshape(sum(PP.a,2),V.h,V.w))
set(gca,'XTick',[],'YTick',[])
colormap gray
ylabel([{'ROI'}],'FontSize',Pl.fs);

% plot F
q=1;
subplot(nrows,ncols,ncols+1), hold on
for j=1:V.Ncells
    F_proj=I{q}.P.a(:,j_inf(j))\F;
    plot(0.9*z1(F_proj(Pl.x_range)),'Color',Pl.c(j,:),'LineWidth',1)
    stem(Pl.n(Pl.x_range,j)-Pl.shift(j),'LineStyle','none','Marker',Pl.m(j),'MarkerEdgeColor',Pl.c(j,:),'MarkerFaceColor',Pl.c(j,:),'MarkerSize',Pl.ms)
    axis([Pl.xlim-Pl.xlim(1) 0 1])
    set(gca,'YTick',[0 1],'YTickLabel',[])
    if q==1
        ylab=ylabel([{'fluorescence'}],'FontSize',Pl.fs);
        set(gca,'XTick',Pl.XTick-min(Pl.XTick),'XTickLabel',(Pl.XTick-min(Pl.XTick))*V.dt,'FontSize',Pl.fs)
        xlabel('time (sec)','FontSize',Pl.fs)
    else
        set(gca,'XTick',Pl.XTick-min(Pl.XTick),'XTickLabel',[])
    end
end

for q=1:2

    % plot spatial filters
    for j=1:V.Ncells,
        subplot(nrows,ncols,1+q+(j-1)*2)
        if q==1, j_inf=[1 2]; else j_inf=[2 1]; end
        imagesc(reshape(EI{q}.P.a(:,j_inf(j)),V.h,V.w)),
        set(gca,'XTickLabel',[],'YTickLabel',[])

        if q==1 && j==1
            ylabel([{'spatial'}; {'filter'}],'FontSize',Pl.fs);
            title('Truth','FontSize',Pl.fs),
        else
            title('Estimated','FontSize',Pl.fs)
        end
    end


    % plot inferred spike train
    for j=1:V.Ncells
        subplot(nrows,ncols,1+q+(j-1)*2+ncols)
        hold on
        if q==1, j_inf=[1 2]; else j_inf=[2 1]; end
        stem(Pl.x_range,Pl.n(Pl.x_range,j)+.03,'LineStyle','none','Marker','v','MarkerEdgeColor',Pl.c(j,:),'MarkerFaceColor',Pl.c(j,:),'MarkerSize',Pl.ms)
        bar(Pl.x_range,I{q}.n(Pl.x_range,j_inf(j))/max(I{q}.n(Pl.x_range,j_inf(j))),'EdgeColor',Pl.c(j,:),'FaceColor',Pl.c(j,:))
        axis('tight')
        set(gca,'YTick',[0 1],'YTickLabel',[])
        set(gca,'XTick',Pl.XTick,'XTickLabel',(Pl.XTick-min(Pl.XTick))*V.dt,'FontSize',Pl.fs)
        xlabel('time (sec)','FontSize',Pl.fs)
        if q==1 && j==1
            ylabel([{'fast'}; {'filter'}],'FontSize',Pl.fs)
        else
        end
    end

end

% neuron={'                 Neuron 1'};
% annotation(gcf,'textbox',[0.2642 0 0.3248 0.6531],'String',neuron,'FitBoxToText','off');

% print fig
wh=[7 2];   %width and height
DirName = '../../figs/';
PrintFig(wh,DirName,fname);
