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

% 1) load data

clear, clc
cd ~/Research/oopsi/fast-oopsi/code/scripts/
load('../../../meta-oopsi/data/tanya/s2m2_40xcorrect.mat')
D=dd;

%%
ind =[
    2
    3
    4
    5
    6
    7
    8
    9
    10
    22
    23
    24
    25
    26
    27
    28
    29
    30
    42
    43
    44
    45
    46
    47
    48
    49
    50
    62
    63
    64
    65
    66
    67
    68
    69
    70
    82
    83
    84
    85
    86
    87
    88
    89
    90
    102
    103
    104
    105
    106
    107
    108
    109
    110
    122
    123
    124
    125
    126
    127
    128
    129
    130
    142
    143
    144
    145
    146
    147
    148
    149
    150
    162
    163
    164
    165
    166
    167
    168
    169
    170
    182
    183
    184
    185
    186
    187
    188
    189
    190
    202
    203
    204
    205
    206
    207
    208
    209
    210
    222
    223
    224
    225
    226
    227
    228
    229
    230
    242
    243
    244
    245
    246
    247
    248
    249
    250
    262
    263
    264
    265
    266
    267
    268
    269
    270
    282
    283
    284
    285
    286
    287
    288
    289
    290];
%%

F=D.F(ind,:)';
n=D.n';
siz=size(F);

% 2) set variables

V.T     = siz(1);                    % # of time steps
V.dt    = 0.0331;                  % time step size
V.Np    = siz(2);                  % # of pixels in each image
V.w     = 9;                  % width of frame (pixels)
V.h     = 15;                 % height of frame (pixels)
V.Ncells= 1;                     % # cells
V.plot  = 0;                      % whether to plot filter with each iteration
V.save  = 0;
V.name  = 'spatial_learn_multi_vitro';

Fmedian=median(F,1)';
F=detrend(F);
F=F-min(F(:));
F=F/max(F(:));

% figure(1), clf, imagesc(reshape(mean(F,1),V.w,V.h)), colorbar


%% 5) make movie

MakMov  = 0;
% make movie of raw data
if MakMov==1
    for i=1:V.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(F(i,:),width,height),['../../data/' V.name],'tif','Compression','none','WriteMode',mod)
    end
end

%% 6) infer spike train using various approaches

P=struct;
qs=[1:3]; %[1 2 3 4];
for q=qs
    GG=F; Tim=V;
    %     if q==1,
    %         Phat{q}=P;
    %         fast{q}.label='true filter';
    if q==1
        Phat{q}=P;
        boxcar=zeros(V.w,V.h);
        boxcar(5:9,7:12)=1;
        Phat{q}.a=boxcar(:);
        Phat{q}.b=mean(F)';
        fast{q}.label='boxcar filter';
    elseif q==2
        Phat{q}=P;
        Phat{q}.a=1+0*boxcar(:); %Fmedian;
        Phat{q}.b=mean(F)';
        fast{q}.label='median filter';
    elseif q==3
        Phat{q}=P;
        Phat{q}.a=Phat{2}.a;
        Phat{q}.b=Phat{2}.b;
        fast{q}.label='learned filter';
        Tim.fast_iter_max=1;
        Tim.fast_ignore_post=1;
        Tim.fast_plot=1;
        Tim.est_a=1;
        Tim.est_b=1;
    elseif q==4
        fast{q}.label='scalarized filter';
        GG=Phat{1}.a\F';
        GG=detrend(GG);
        GG=GG-min(GG);
        GG=GG/max(GG);
        Tim.Np=1;
        Tim.Ncells=1;
        Phat{q}=struct;
        Phat{q}.lam=1;
    end
    display(fast{q}.label)
    starttime=cputime;
    [fast{q}.n fast{q}.P fast{q}.V] = fast_oopsi(GG',Tim,Phat{q});
    fast{q}.V.time = cputime-starttime;
end
if V.save==1, save(['../../data/' V.name]); end

%% plot results
% V.name  = 'spatial_inf';
% load(['../../data/' V.name])

nrows   = 3;                                 % set number of rows
ncols   = length(qs);
h       = zeros(nrows,1);
xlims   = [1 V.T];                            % time steps to plot
xs      = xlims(1):xlims(2);
nticks  = 5;                                    % number of ticks along x-axis
fs      = 14;                       % font size
ms      = 5;                        % marker size for real spike
sw      = 2;                        % spike width
lw      = 2;                        % line width
Nc      = V.Ncells;
XTicks  = [xlims(1):range(xlims)/4:xlims(2)]-xlims(1);
gray    = 0.4*[1 1 1];            % define gray color

h       = zeros(nrows*ncols);
fig     = figure(1); clf,

ntemp=double(n);
ntemp(ntemp==0)=NaN;


for q=qs
    if q==1, i=0; elseif q==2 i=1; elseif q==3, i=2; 
    elseif q==4, i=3; Phat{q}=Phat{1}; end
    
    % plot filter
    i=i+1;
    subplot(nrows,ncols,i)
    imagesc(reshape(Phat{q}.a,V.w,V.h))
    colormap('gray')
    title(fast{q}.label,'FontSize',fs)
    if q==1,
        ylab=ylabel([{'spatial'}; {'filters'}],'FontSize',fs);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    end
    
    % plot F projection
    i=i+ncols;
    h(i)=subplot(nrows,ncols,i);
    plot(z1(Phat{q}.a\F(xs,:)'),'k','LineWidth',lw)
    hold on
    stem(ntemp(xs)+0.1,'Marker','+','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)

    set(gca,'XTick',XTicks,'XTickLabel',[],'FontSize',fs)
    set(gca,'YTick',[],'FontSize',fs)
    axis([xlims-xlims(1) 0 1.1])
    if q==1,
        ylab=ylabel([{'fluorescence'}],'FontSize',fs);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    end
    box off

    % plot nhat
    i=i+ncols;
    h(i)=subplot(nrows,ncols,i);
    bar(fast{q}.n(xs))
    hold on
    stem(z1(ntemp(xs)),'Marker','+','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
    axis('tight')
    set(gca,'XTick',XTicks,'XTickLabel',round(XTicks*V.dt),'FontSize',fs)
    set(gca,'YTick',[],'FontSize',fs)
    xlabel('time (sec)','FontSize',fs)
    if q==1,
        ylab=ylabel([{'fast'};  {'filter'}],'FontSize',fs);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    end
    box off
end

ind=find(h==0);
h(ind)=[];
linkaxes(h,'x')

if V.save==1 % print fig
    wh=[10 5];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=['../../figs/' V.name];
    print('-depsc',figname)
    print('-dpdf',figname)
    saveas(fig,figname)
end