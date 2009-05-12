% this script is a wrapper that loads some data (or simulates it) and then
% calls the most recent FOOPSI function to infer spike trains and
% parameters

clear, clc, 
fname = 'spatial_multi';

%% set parameters 

% generate spatial filters
Nc      = 1;                                % # of cells in the ROI
neur_w  = 10;                               % height per neuron
height  = 15;                               % height of frame (pixels)
width   = Nc*neur_w;                        % width of frame (pixels)
Npixs   = height*width;                     % # pixels in ROI
x       = linspace(-5,5,width);
y       = linspace(-5,5,height);
[X,Y]   = meshgrid(x,y);
g1      = zeros(Npixs,Nc);
g2      = 0*g1;
Sigma1  = diag([1,1])*3;                    % var of positive gaussian
Sigma2  = diag([1,1])*5;                    % var of negative gaussian
mu      = [1 1]'*linspace(-1.4,1.4,Nc);     % means of gaussians for each cell (distributed across pixel space)
w       = ones(Nc,1); %Nc:-1:1;                          % weights of each filter
for i=1:Nc
    g1(:,i)  = w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma1);
    g2(:,i)  = 0*w(i)*mvnpdf([X(:) Y(:)],mu(:,i)',Sigma2);
end

% set simulation metadata
Meta.T       = 2000;                        % # of time steps
Meta.dt      = 0.005;                       % time step size
Meta.Np      = Npixs;                       % # of pixels in each image
Meta.h       = height;                      % height of frame (pixels)
Meta.w       = width;                       % width of frame (pixels)

% initialize params
for i=1:Nc
    P.a(:,i)=g1(:,i)-g2(:,i);
end
P.b     = 0.1*ones(Nc,1);                              % baseline is a scaled down version of the sum of the spatial filters

P.sig   = 0.01;                             % stan dev of noise (indep for each pixel)
C_0     = 0;                                % initial calcium
tau     = round(100*rand(Nc,1))/100+0.05;   % decay time constant for each cell
P.gam   = 1-Meta.dt./tau(1:Nc);             % set gam
P.lam   = linspace(15,30,Nc)'.*ones(Nc,1);                     % rate-ish, ie, lam*dt=# spikes per second

%% simulate data
n=zeros(Meta.T,Nc);                         % pre-allocate memory for spike train
C=n;                                        % pre-allocate memory for calcium
for i=1:Nc
    n(1,i)      = C_0;
    n(2:end,i)  = poissrnd(P.lam(i)*Meta.dt*ones(Meta.T-1,1));  % simulate spike train
    n(n>1)      = 1;                                            % make only 1 spike per bin
    C(:,i)      = filter(1,[1 -P.gam(i)],n(:,i));               % calcium concentration
end
Z = 0*n(:,1);
F = (C+repmat(((1-P.gam).*P.b)',Meta.T,1))*P.a' + P.sig*randn(Meta.T,Npixs);              % fluorescence

% set user defined parameters
User.MaxIter = 25;                          % # iterations of EM to estimate params
User.Nc      = Nc;                          % # cells per ROI
User.Plot    = 1;                           % whether to plot filter with each iteration
User.Thresh  = 1;                           % whether to threshold spike train before estimating params (we always keep this on)
User.n       = n;                           % true spike train (only included if known)

save([fname '.mat'],'F','n','P','Meta','User')

%% infer spike trains and parameters

for q=1:2
    if q==1;                                        % use true params
        PP=P;
        User.MaxIter = 1;
    elseif q==2                                     % use svd spatial filter
        [U,S,V]=pca_approx(F,User.Nc);
        PP=P;
        for j=1:User.Nc, PP.a(:,j)=V(:,j); end
%         User.MaxIter = 0;
%     elseif q==3                                     % estimate params
%         PP=P;
%         for j=1:User.Nc, PP.a(:,j)=V(:,j); end
        User.MaxIter = 20;
    end
    [I{q}.n I{q}.P] = FOOPSI_v3_03_01(F,PP,Meta,User);
end

save([fname '.mat'],'-append','I')

%% plot results

Pl.g    = [0.75 0.75 0.75];     % gray color
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
Pl.n    = n; Pl.n(Pl.n==0)=NaN; % true spike train (0's are NaN's so they don't plot)
Pl.T    = Meta.T;
Pl.c    = [0 0 0; Pl.g; 1 1 1]; % colors: black, grey, white
Pl.xlim = [200 300]+14;        % limits of x-axis
Pl.x_range = Pl.xlim(1):Pl.xlim(2);
Pl.XTick= [Pl.xlim(1) round(mean(Pl.xlim)) Pl.xlim(2)];
Pl.XTickLabel = round((Pl.XTick-min(Pl.XTick))*Meta.dt*100)/100;

%% show how our estimation procedure given truth and when estimating spatial filter

figure(1), clf, hold on
nrows   = 2*length(I);
ncols   = 1+User.Nc;

Fmean=mean(F);
for q=1:length(I)
    E{q}.a=I{q}.P.a;
    E{q}.b=I{q}.P.b;
end

maxx=[]; minn=[];
for q=1:length(I)
    E{q}.a_max=max(E{q}.a(:));
    E{q}.a_min=min(E{q}.a(:));
    for j=1:User.Nc
        E{q}.a(:,j)=60*(E{q}.a(:,j)-E{q}.a_min)/(E{q}.a_max-E{q}.a_min);
    end
end

for q=1:length(I)

    % align inferred cell with actual one
    j_inf=0*n(1:User.Nc);
    for j=1:User.Nc
        cc=0*n(1:User.Nc);
        for k=1:User.Nc
            cc_temp=corrcoef(n(:,j),I{q}.n(:,k)); cc(k)=cc_temp(2);
        end
        [foo j_inf(j)]=max(cc);
    end

    % plot mean image frame
    subplot(nrows,ncols, (q-1)*2*ncols+1)
    imagesc(reshape(sum(E{q}.a,2),Meta.h,Meta.w))
    set(gca,'XTick',[],'YTick',[])
    colormap gray
    if q==1, title(['Sum of Spatial Filters'],'FontSize',Pl.fs); end

    % plot spatial filters
    for j=1:Nc,
        subplot(nrows,ncols,(q-1)*2*ncols+j+1)
        image(reshape(E{q}.a(:,j),Meta.h,Meta.w)),
        set(gca,'XTickLabel',[],'YTickLabel',[])
        if q==1, title(['Neuron ' num2str(j)],'FontSize',Pl.fs); end
    end

    % plot F and n for neurons
    subplot(nrows,ncols,(q-1)*2*ncols+ncols+1), hold on
    for j=1:User.Nc
        F_proj=E{q}.a(:,j)\(F)';
        plot(z1(F_proj(Pl.x_range))+1,'Color',Pl.c(j,:),'LineWidth',Pl.lw)
        bar(Pl.n(Pl.x_range,j),'FaceColor',Pl.c(j,:),'EdgeColor',Pl.c(j,:),'LineWidth',2)
        axis([Pl.xlim-Pl.xlim(1) 0 2])
        set(gca,'XTick',Pl.XTick,'XTickLabel',[])
        set(gca,'YTick',[0 1 2.5],'YTickLabel',[])
        set(gca,'XTick',Pl.XTick,'XTickLabel',[])
    end

    % plot inferred spike train
    for j=1:User.Nc
        subplot(nrows,ncols,(q-1)*2*ncols+ncols+1+j)
        hold on
        stem(Pl.x_range,Pl.n(Pl.x_range,j),'LineStyle','none','Marker','v','MarkerEdgeColor',Pl.g,'MarkerFaceColor',Pl.g,'MarkerSize',Pl.ms)
        if j==1, k=2; else k=1; end
        bar(Pl.x_range,I{q}.n(Pl.x_range,j_inf(j))/max(I{q}.n(Pl.x_range,j_inf(j))))
        axis('tight')
        set(gca,'YTick',[0 1],'YTickLabel',[])
        set(gca,'XTick',Pl.XTick-min(Pl.XTick),'XTickLabel',[])
    end

end

% print fig
wh=[7.5 7];   %width and height
DirName = '../../docs/journal_paper/figs/';
FileName = 'multi_spacial';
PrintFig(wh,DirName,FileName);

%% make hist stuff

edges=linspace(0,1,20);
for q=1:length(I)
    D{q}.n      = I{q}.n;
    D{q}.sp     = [];
    D{q}.no_sp  = [];
    for j=1:User.Nc
        D{q}.n(:,j)=I{q}.n(:,j)/max(I{q}.n(:,j));

        D{q,j}.sp       = D{q}.n(:,j);
        D{q,j}.no_sp    = D{q,j}.sp; %D{q}.n(:,j);

        D{q,j}.sp(n(:,j_inf(j))==0)=[];
        D{q}.sp=[D{q}.sp; D{q,j}.sp];


        D{q,j}.no_sp(n(:,j_inf(j))==1)=[];
        D{q}.no_sp=[D{q}.no_sp; D{q,j}.no_sp];

        D{q,j}.hist(:,1)=histc(D{q,j}.no_sp,edges);
        D{q,j}.hist(:,2)=histc(D{q,j}.sp,edges);
    end
end


%% multi_hist1
figure(9), clf, hold on

ncols=2;%length(I);
nrows=2;
xticks=[0 .5 1];
for q=1:length(I)
    for k=1:2
        joint_hist=[];

        for j=1:User.Nc
            joint_hist=[joint_hist D{q,j}.hist(:,k)];
        end
        joint_hist=joint_hist./repmat(sum(joint_hist),length(joint_hist),1);
        
        subplot(nrows,ncols,q+(k-1)*nrows)
        hold on

        for j=1:User.Nc
            bar(edges+diff(edges(1:2)/User.Nc),joint_hist(:,j),'FaceColor',Pl.c(j,:),'EdgeColor',Pl.c(j,:),'BarWidth',1/User.Nc)
        end
%         bar(edges,joint_hist(:,1),'FaceColor','k','EdgeColor','k','BarWidth',0.5)

        set(gca,'XTick',xticks,'XTickLabel',xticks,'YTickLabel',[])
        axis('tight')
        if q==1 && k==1
            ylab=ylabel(['P(n_{MAP} | n=1)'],'FontSize',Pl.fs);
            title('Optimal Spatial Filters','FontSize',Pl.fs)
        elseif q==1 && k==2
            ylab=ylabel(['P(n_{MAP} | n=0)'],'FontSize',Pl.fs);

        elseif q==2 && k==1
            title('Estimated Spatial Filters','FontSize',Pl.fs)
        end
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    end
end

% print fig
wh=[7.5 7];   %width and height
DirName = '../../docs/journal_paper/figs/';
FileName = 'multi_hist1';
PrintFig(wh,DirName,FileName);

%% plot distribution of spikes 

% edges=linspace(0,1,20);
% for q=1:length(I)
%     D{q}.n=I{q}.n;
%     D{q}.sp=[];
%     D{q}.no_sp=[];
%     for j=1:User.Nc
%         D{q}.n(:,j)=D{q}.n(:,j)/max(D{q}.n(:,j));
% 
%         D{q,j}.sp=D{q}.n(:,j);
%         D{q,j}.sp(n(:,j_inf(j))==0)=[];
%         D{q}.sp=[D{q}.sp; D{q,j}.sp];
% 
% 
%         D{q,j}.no_sp=D{q}.n(:,j);
%         D{q,j}.no_sp(n(:,j_inf(j))==1)=[];
%         D{q}.no_sp=[D{q}.no_sp; D{q,j}.no_sp];
%     end
% end
% 
% xspace=.207;
% yspace=.285;
% w=xspace/2;
% h=yspace/2;
% 
% for q=1:length(I)
%     for j=1:User.Nc
%         subplot('Position',[0.75 0.68+h .157/2 h]), 
%         D{q,j}.hist=histc(D{q,j}.no_sp,edges);
%         bar(edges,D{q,j}.hist)
%         set(gca,'XTick',[0 .5 1],'XTickLabel',[],'YTickLabel',[])
% 
%         %    subplot(nrows,ncols,7+0.4*j)
%         %    D{q,j}.hist=histc(D{q,j}.sp,edges);
%         %    bar(edges,D{q,j}.hist)
%         %    set(gca,'XTickLabel',[0 .5 1],'XTickLabel',[],'YTickLabel',[])
%     end
% end
% 
% %%
% figure(9),
% for q=1:length(I)
%     D{q}.sp=[];
%     D{q}.no_sp=[];
%     for j=1:User.Nc
% 
%         D{q,j}.sp=D{q}.n(:,j);
%         D{q,j}.sp(n(:,j_inf(j))==0)=NaN;
%         D{q}.sp=[D{q}.sp; D{q,j}.sp];
% 
% 
%         D{q,j}.no_sp=D{q}.n(:,j);
%         D{q,j}.no_sp(n(:,j_inf(j))==1)=NaN;
%         D{q}.no_sp=[D{q}.no_sp; D{q,j}.no_sp];
%     end
% end
% boxplot([D{q}.sp D{q}.no_sp D{2}.sp D{2}.no_sp D{3}.sp D{3}.no_sp],'notch','on','whisker',115)
