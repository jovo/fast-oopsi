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
F=D.F';
n=D.n';
siz=size(F);

% 2) set variables

V.T     = siz(1);                    % # of time steps
V.dt    = 0.0331;                  % time step size
V.Np    = siz(2);                  % # of pixels in each image
V.w     = D.nRows;                  % width of frame (pixels)
V.h     = D.nCols;                 % height of frame (pixels)
V.Ncells= 2;                     % # cells
V.plot  = 0;                      % whether to plot filter with each iteration
V.save  = 1;
V.name  = 'spatial_learn_multi_vitro';
V.fast_iter_max=1;

V.T=1000;
F0=F;
F=F(1:V.T,:);
F=detrend(F')';
F=F-min(F(:));
F=F/max(F(:));


figure(1), clf, imagesc(reshape(mean(F,1),V.w,V.h)), colorbar

%% 3) infer spike train using various approaches
clf
qs=[1];

boxcar{1}=zeros(V.w,V.h);
boxcar{1}(6:11,6:10)=1;
boxcar{2}=zeros(V.w,V.h);
boxcar{2}(14:18,12:17)=1;
Phat{1}.a=[boxcar{1}(:) boxcar{2}(:)];
Phat{1}.b=0*boxcar{1}(:);

GG=F;
GG=GG-repmat(mean(GG),V.T,1);
  
for q=qs
    GGG=0*GG;
            
%     GG=GG./repmat(max(GG),V.T,1);
%     GG=detrend(GG')';
    Tim=V;
    if q==1,
        fast{q}.label='boxcar filter';
        display(fast{q}.label)
        starttime=cputime;
        
        for i=1:V.Np
            f = detrend(GG(:,i));
            nfft    = 2^nextpow2(V.T);
            y       = fft(f,nfft);
            y(round(295/2+1):round(305/2+1)) = 0;
            bw      = 10;
            y(1:bw) = 0; y(end-bw+1:end)=0;
            iy      = ifft(y,nfft);
            GGG(:,i)= z1(real(iy(1:V.T)));
        end
        GGG=GGG-repmat(mean(GGG),V.T,1);


        [fast{q}.n fast{q}.P fast{q}.V] = fast_oopsi(GGG',Tim,Phat{q});
        fast{q}.V.time = cputime-starttime;
    elseif q==2
        fast{q}.label='scalarized filter';
        GGG=Phat{1}.a\F';
        Tim.Np=1;
        Tim.Ncells=1;
        Tim.est_sig=0;
        Tim.fast_iter_max=1;
        Tim.fast_plot=1;
        Phat{q}=struct;
        Phat{q}.lam=1;
        display(fast{q}.label)
        starttime=cputime;
        for i=1:V.Ncells

            f = detrend(GGG(i,:));
            nfft    = 2^nextpow2(V.T);
            y       = fft(f,nfft);
            y(round(295/2+1):round(305/2+1)) = 0;
            bw      = 10;
            y(1:bw) = 0; y(end-bw+1:end)=0;
            iy      = ifft(y,nfft);
            GGG(i,:)= z1(real(iy(1:V.T)));
            
            [ftemp{i}.n ftemp{i}.P ftemp{i}.V] = fast_oopsi(GGG(i,:),Tim,Phat{q});
        end
        fast{q}.n = [ftemp{1}.n ftemp{2}.n];
        fast{q}.P.a = Phat{1}.a;
        fast{q}.V.time = cputime-starttime;
    elseif q==3

        for i=1:V.Np
            f = detrend(GG(:,i));
            nfft    = 2^nextpow2(V.T);
            y       = fft(f,nfft);
            y(round(295/2+1):round(305/2+1)) = 0;
            bw      = 10;
            y(1:bw) = 0; y(end-bw+1:end)=0;
            iy      = ifft(y,nfft);
            GGG(:,i)= z1(real(iy(1:V.T)));
        end
        GGG=GGG-repmat(mean(GGG),V.T,1);
        
        Phat{q}=Phat{1};
        fast{q}.label='learned filter';
        Tim.fast_iter_max=15;
        Tim.fast_ignore_post=1;
        Tim.fast_plot=1;
        Tim.est_a=1;
        Tim.est_b=1;
        Tim.est_sig=0;
        Tim.est_lam=0;
        Tim.est_gam=0;
        [fast{q}.n fast{q}.P fast{q}.V] = fast_oopsi(GGG',Tim,Phat{q});

    elseif q==4
        Phat{q}=Phat{1};
%         Phat{q}.b=-mean(F)';
        fast{q}.label='learned filter';
        Tim.fast_iter_max=15;
        Tim.fast_ignore_post=1;
        Tim.fast_plot=1;
        Tim.est_a=1;
        Tim.est_b=1;
        Tim.est_sig=0;
        Tim.est_lam=0;
        Tim.est_gam=0;
        [fast{q}.n fast{q}.P fast{q}.V] = fast_oopsi(GG',Tim,Phat{q});
    elseif q==6
        Phat{q}=Phat{1};
%         Phat{q}.b=-mean(F)';
        fast{q}.label='learned filter';
        Tim.fast_iter_max=15;
        Tim.fast_ignore_post=1;
        Tim.fast_plot=1;
        Tim.est_a=1;
        Tim.est_b=0;
        Tim.est_sig=0;
        Tim.est_lam=0;
        Tim.est_gam=0;
        [fast{q}.n fast{q}.P fast{q}.V] = fast_oopsi(GG',Tim,Phat{q});
    elseif q==7
                for i=1:V.Np
            f = detrend(GG(:,i));
            nfft    = 2^nextpow2(V.T);
            y       = fft(f,nfft);
            y(round(295/2+1):round(305/2+1)) = 0;
            bw      = 10;
            y(1:bw) = 0; y(end-bw+1:end)=0;
            iy      = ifft(y,nfft);
            GGG(:,i)= z1(real(iy(1:V.T)));
        end
        GGG=GGG-repmat(mean(GGG),V.T,1);
        
        Phat{q}=Phat{1};
        fast{q}.label='learned filter';
        Tim.fast_iter_max=15;
        Tim.fast_ignore_post=1;
        Tim.fast_plot=1;
        Tim.est_a=1;
        Tim.est_b=0;
        Tim.est_sig=0;
        Tim.est_lam=0;
        Tim.est_gam=0;
        [fast{q}.n fast{q}.P fast{q}.V] = fast_oopsi(GGG',Tim,Phat{q});
    
    end
end
if V.save==1, save(['../../data/' V.name]); end

%% plot results
% V.name  = 'spatial_learn_multi_vitro';
% load(['../../data/' V.name])

qs=[1 3 7];
nrows   = 6;                                 % set number of rows
ncols   = length(qs);
h       = zeros(nrows,1);
xlen    = 270;
xlims   = [0 500]; %+300; %round([V.T-xlen V.T+xlen]/2);                            % time steps to plot
xs      = xlims(1)+1:xlims(2);
nticks  = 4;                                    % number of ticks along x-axis
fs      = 14;                       % font size
ms      = 5;                        % marker size for real spike
sw      = 2;                        % spike width
lw      = 2;                        % line width
Nc      = V.Ncells;
XTicks  = [xlims(1):180:xlims(2)]-xlims(1);
gray    = 0.7*[1 1 1];            % define gray color

fig = figure(2); clf,


subplot(nrows,ncols,1)
imagesc(reshape(sum(-F0,1),V.w,V.h))
title('mean frame','FontSize',fs)
ylab=ylabel([{'summed'}; {'spatial'}; {'filters'}],'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')


subplot(nrows,ncols,2)
imagesc(reshape(-F0(round(mean(xlims)),:),V.w,V.h))
title('typical frame','FontSize',fs)

subplot(nrows,ncols,3)
imagesc(reshape(fast{3}.P.b,V.w,V.h))
title('estimated b','FontSize',fs)

ntemp=double(n(xs,:));
ntemp(ntemp==0)=NaN;
fast{1}.label='boxcar'; fast{4}.label='learned';

for q=qs
    if q==1, i=0; elseif q==3, i=1; elseif q==7, i=2; end    
%     i=q-1;
    i=i+ncols;
    
    i=i+1;
    subplot(nrows,ncols,i)
    imagesc(reshape(sum(fast{q}.P.a,2),V.w,V.h))
    colormap('gray')
    title(fast{q}.label,'FontSize',fs)
    
    i=i+ncols;
    subplot(nrows,ncols,i)
    imagesc(reshape(fast{q}.P.a(:,1),V.w,V.h))
    colormap('gray')
    if q==1,
        ylab=ylabel([{'spatial'}; {'filter'}],'FontSize',fs);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    end

    i=i+ncols;
    subplot(nrows,ncols,i)
    FF=(fast{q}.P.a\F(xs,:)')';
    hold on
    plot(z1(FF(:,1))+1.2,'k','LineWidth',lw)
    stem(ntemp(:,1)+0.1,'Marker','+','MarkerSize',ms,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
    bar(z1(fast{q}.n(xs,1)))
    
    set(gca,'XTick',XTicks,'XTickLabel',[],'FontSize',fs)
    set(gca,'YTick',[],'FontSize',fs)
    axis([xlims-xlims(1) 0 2.2])

    if q==1,
        ylab=ylabel([{'fluorescence'}; {'fast filter'}],'FontSize',fs);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    end
    box off

    i=i+ncols;
    subplot(nrows,ncols,i)
    imagesc(reshape(fast{q}.P.a(:,2),V.w,V.h))
    colormap('gray')
    if q==1,
        ylab=ylabel([{'spatial'}; {'filter'}],'FontSize',fs);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    end

    i=i+ncols;
    subplot(nrows,ncols,i)
    hold on

    plot(z1(FF(:,2))+1.2,'color','k','LineWidth',lw)
    stem(ntemp(:,2)+0.1,'Marker','+','MarkerSize',ms,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
    bar(z1(fast{q}.n(xs,2)),'EdgeColor','k','FaceColor','k')

    axis([xlims-xlims(1) 0 2.2])
    set(gca,'XTick',XTicks,'XTickLabel',round(XTicks*V.dt),'FontSize',fs)
    set(gca,'YTick',[],'FontSize',fs)
    xlabel('time (sec)','FontSize',fs)
    if q==1,
        ylab=ylabel([{'fluorescence'}; {'fast filter'}],'FontSize',fs);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    end
    box off
end


l1=.79; %0.7631;
l2=.65; %0.4236;
l3=.35;

% annotation(gcf,'line',[0 0.91],[l1 l1]);
annotation(gcf,'line',[0 0.91],[l2 l2]);
annotation(gcf,'line',[0 0.91],[l3 l3]);

annotation(gcf,'textbox',[0.0 0.9348 0.04713 0.03339],'String','(a)','FitBoxToText','on','LineStyle','none','FontSize',fs);
% annotation(gcf,'textbox',[0.0 0.7600 0.04713 0.03339],'String','(b)','FitBoxToText','on','LineStyle','none','FontSize',fs);
annotation(gcf,'textbox',[0.0 0.6200 0.04713 0.03339],'String','(b)','FitBoxToText','on','LineStyle','none','FontSize',fs);
annotation(gcf,'textbox',[0.0 0.3200 0.04713 0.03339],'String','(c)','FitBoxToText','on','LineStyle','none','FontSize',fs);
            
if V.save==1 % print fig
    wh=[9 6];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=['../../figs/' V.name];
    print('-depsc',figname)
    print('-dpdf',figname)
    saveas(fig,figname)
end