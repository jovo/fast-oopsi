% this script is a wrapper that loads some data (or simulates it) and then
% calls the most recent FOOPSI function to infer spike trains and
% parameters

clear, clc,
fname = 'spatial_multi_learn';

%% set parameters

% generate spatial filters
Nc      = 2;                                % # of cells in the ROI
neur_w  = 8;                                % height per neuron
height  = 10;                               % height of frame (pixels)
width   = Nc*neur_w;                        % width of frame (pixels)
Npixs   = height*width;                     % # pixels in ROI
x       = linspace(-5,5,width);
y       = linspace(-5,5,height);
[X,Y]   = meshgrid(x,y);
g1      = zeros(Npixs,Nc);
g2      = 0*g1;
Sigma1  = diag([1,1])*3;                    % var of positive gaussian
Sigma2  = diag([1,1])*5;                    % var of negative gaussian
sp      = 1.8;
mu      = [1 1]'*linspace(-sp,sp,Nc);       % means of gaussians for each cell (distributed across pixel space)
w       = linspace(1,Nc,Nc); %Nc:-1:1;      % weights of each filter
for i=1:Nc
    g1(:,i)  = w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma1);
    g2(:,i)  = 0*w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma2);
end

for i=1:Nc
    P.a(:,i)=g1(:,i)-g2(:,i);
end
P.b     = 0.01*mean(P.a);               % baseline is a scaled down version of the sum of the spatial filters

% set Variables 
V.T      = 4800;                        % # of time steps
V.dt     = 1/60;                        % time step size
V.Npixels= Npixs;                       % # of pixels in each image
V.Ncells = Nc;                          % # cells per ROI
V.h       = height;                     % height of frame (pixels)
V.w       = width;                      % width of frame (pixels)

% initialize parameters for estimating parameters (these are only used if MaxIter>1)
V.est_sig     = 0;
V.est_lam     = 0;
V.est_gam     = 0;
V.est_b       = 1;
V.est_a       = 1;

% initialize parameters for plotting results after each pseudo-EM step
V.fast_thresh  = 1;                     % whether to threshold spike train before estimating 'a' and 'b' (we always keep this on)
V.fast_plot    = 1;                     % whether to plot filter with each iteration
V.fast_poisson = 0;                     % whether observations are poisson or gaussian
V.fast_iter_max= 3;                     % # iterations of EM to estimate params
V.fast_ignore_post=1;


% initialize other parameters
P.sig   = max(P.a(:))*1;               % stan dev of noise (indep for each pixel)
C_0     = 0;                             % initial calcium
tau     = rand(V.Ncells,1)/2+.05;        % decay time constant for each cell
P.gam   = 1-V.dt./tau(1:Nc);             % set gam
lam     = linspace(0,4*pi,V.T);
P.lam = repmat(exp(1.6*(sin(lam))'),1,V.Ncells); % rate-ish, ie, lam*dt=# spikes per second
for j=2:V.Ncells
    P.lam(:,j) = exp(1.6*(sin(lam+pi))');
end
P.lam   = P.lam-repmat(min(P.lam),V.T,1);

%% simulate data

V.n=zeros(V.T,V.Ncells);                   % pre-allocate memory for spike train
V.C=V.n;                                     % pre-allocate memory for calcium
for i=1:V.Ncells
    V.n(1,i)      = C_0;
    V.n(2:end,i)  = poissrnd(P.lam(2:end,i)*V.dt);    % simulate spike train
    V.n(V.n>1)      = 1;                                % make only 1 spike per bin
    V.C(:,i)      = filter(1,[1 -P.gam(i)],V.n(:,i));   % calcium concentration
end
V.F = P.a*(V.C+repmat(P.b,V.T,1))'+P.sig*rand(Npixs,V.T);
V.n(V.n==0)=NaN;
save(['../../data/' fname '.mat'],'P','V')

%% generate tif

MakMov  = 1;
if MakMov==1
    FF=uint8(floor(255*z1(V.F)));
    for t=1:500
        if t==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(FF(:,t),height,width),['../../data/' fname '.tif'],'tif','Compression','none','WriteMode',mod)
    end
end

%% infer spike trains and parameters
F=V.F-repmat(mean(V.F,2),1,V.T); 
n=V.n;
V.est_b=1;
V.fast_ignore_post=1;
V.fast_iter_max=5;

% infer spike trains using a variety of techniques
q=1;
exps=1;
PP=P;
PP.lam=mean(PP.lam);
[U,S,VV]=pca_approx(F',V.Ncells);
for j=1:V.Ncells, PP.a(:,j)=VV(:,j); end
PP.b=0*P.b;%1e-3*ones(1,V.Ncells);
for q=1:2
if q==2, V.fast_thr=1; end
[I{q}.n I{q}.P] = fast_oopsi(F,V,PP);
end
save(['../../data/' fname '.mat'],'-append','I','PP')
sound(10*sin(linspace(0,180*pi,2000)))

%% plot results
% fname = 'spatial_multi_learn';
% load(['../../data/' fname])

Pl.g    = 0.65*ones(1,3);       % gray color
Pl.fs   = 12;                   % font size
Pl.w1   = 0.28;                 % width of subplot
Pl.wh   = [Pl.w1 Pl.w1];        % width and height of subplots
Pl.b1   = 0.67;                 % bottom of topmost subplt
Pl.l1   = 0.41;                 % left side of left subplot on right half
Pl.l2   = Pl.l1+0.3;            % left side of right subplot on right half
Pl.s1   = Pl.w1+0.043;          % space between subplots on right side
Pl.s2   = .33;                  % space between subplots on left side
Pl.ms   = 5;                    % marker size
Pl.lw   = 2;                    % line width
Pl.n    = V.n; Pl.n(Pl.n==0)=NaN; % true spike train (0's are NaN's so they don't plot)
Pl.T    = V.T;
Pl.c    = [0 0 0; Pl.g; 1 1 1]; % colors: black, grey, white
Pl.m    = ['+','+'];
Pl.xlim = [200 700]+2000;        % limits of x-axis
Pl.shift= [0 .07];
Pl.x_range = Pl.xlim(1):Pl.xlim(2);
Pl.XTick   = (linspace(Pl.xlim(1),Pl.xlim(2),5));
Pl.XTickLabel = ceil((Pl.XTick-min(Pl.XTick))*V.dt*100)/100;

% show how our estimation procedure given truth and when estimating spatial filter

figure(1), clf, hold on
nrows   = 2; %*length(I);
ncols   = 1+V.Ncells;

Fmean=mean(F);
for q=1:length(exps)
    E{q}.a=I{q}.P.a;
    E{q}.b=I{q}.P.b;
end

maxx=[]; minn=[];
for q=1:length(exps)
    E{q}.a_max=max(E{q}.a(:));
    E{q}.a_min=min(E{q}.a(:));
    for j=1:V.Ncells
        E{q}.a(:,j)=60*(E{q}.a(:,j)-E{q}.a_min)/(E{q}.a_max-E{q}.a_min);
    end
end

for q=1%:length(exps)

    % align inferred cell with actual one
    j_inf=0*n(1:V.Ncells);
    cc=0*n(1:V.Ncells,:);
    for j=1:V.Ncells
        for k=1:V.Ncells
            cc_temp = corrcoef(n(:,j),I{q}.n(:,k));
            cc(k,j)   = cc_temp(2);
        end
    end
    [foo ind] = max(cc);
    sortI = sort(ind);
    if ~any(diff(sortI)==0)
        j_inf=ind;
    else
        [foo ind]=max(cc(:));
        j_inf(1)=(ind+1)/2;
        if j_inf(1)==1; j_inf(2)=2;
        else j_inf(1)=2; j_inf(2)=1; end
    end

    % plot mean image frame
    subplot(nrows,ncols, (q-1)*2*ncols+1)
    imagesc(reshape(sum(V.F,2),V.h,V.w))
    set(gca,'XTick',[],'YTick',[])
    colormap gray
    if q==1,
        title(['mean frame'],'FontSize',Pl.fs);
%         ylab=ylabel([{'truth'}]);
    elseif q==2
        ylab=ylabel([{'PCA'}]);
    elseif q==3
        ylab=ylabel([{'Estimated'}]);
    end
%     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',Pl.fs)

    % plot spatial filters
    for j=1:V.Ncells,
        subplot(nrows,ncols,(q-1)*2*ncols+j+1)
        image(reshape(E{q}.a(:,j_inf(j)),V.h,V.w)),
        set(gca,'XTickLabel',[],'YTickLabel',[])
        if q==1, title(['neuron ' num2str(j)],'FontSize',Pl.fs); end
        ylabel([{'estimated'}; {'spatial filter'}],'FontSize',Pl.fs)
    end

    % plot F and n for neurons
    subplot(nrows,ncols,(q-1)*2*ncols+ncols+1), hold on
    for j=1:V.Ncells
        F_proj=E{q}.a(:,j_inf(j))\F;
        plot(0.9*z1(F_proj(Pl.x_range)),'Color',Pl.c(j,:),'LineWidth',Pl.lw)
        stem(Pl.n(Pl.x_range,j)-Pl.shift(j),'LineStyle','none','Marker',Pl.m(j),'MarkerEdgeColor',Pl.c(j,:),'MarkerFaceColor',Pl.c(j,:),'MarkerSize',Pl.ms)
        axis([Pl.xlim-Pl.xlim(1) 0 1])
        set(gca,'YTick',[0 1],'YTickLabel',[])
        if q==1
            ylabel('fluorescence','FontSize',Pl.fs)
            set(gca,'XTick',Pl.XTick-min(Pl.XTick),'XTickLabel',round((Pl.XTick-min(Pl.XTick))*V.dt))
            xlabel('time (sec)','FontSize',Pl.fs)
        else
            set(gca,'XTick',Pl.XTick-min(Pl.XTick),'XTickLabel',[])
        end
    end

    % plot inferred spike train
    for j=1:V.Ncells
        subplot(nrows,ncols,(q-1)*2*ncols+ncols+1+j)
        hold on
        stem(Pl.x_range,Pl.n(Pl.x_range,j)+.03,'LineStyle','none','Marker','+','MarkerEdgeColor',Pl.c(j,:),'MarkerFaceColor',Pl.c(j,:),'MarkerSize',Pl.ms)
        if j==1, k=2; else k=1; end
        bar(Pl.x_range,I{q}.n(Pl.x_range,j_inf(j))/max(I{q}.n(Pl.x_range,j_inf(j))),'EdgeColor',Pl.c(j,:),'FaceColor',Pl.c(j,:))
        axis('tight')
        set(gca,'YTick',[0 1],'YTickLabel',[])
        if q==1
            ylabel('fast filter','FontSize',Pl.fs)
            set(gca,'XTick',Pl.XTick,'XTickLabel',round((Pl.XTick-min(Pl.XTick))*V.dt))
            xlabel('time (sec)','FontSize',Pl.fs)
        else
            set(gca,'XTick',Pl.XTick,'XTickLabel',[])
        end
    end

end

% print fig
wh=[7.5 3.5*length(exps)];   %width and height
DirName = '../../figs/';
PrintFig(wh,DirName,fname);
