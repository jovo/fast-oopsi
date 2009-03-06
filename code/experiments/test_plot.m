clear Pl
nrows   = 3;                                  % set number of rows
ncols   = 1+Nc;
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
Pl.XTicks=[200 400 600];

fnum=figure(3); clf

i=1; h(i) = subplot(nrows,ncols,i);
imagesc(reshape(sum(P.a'),Sim.w,Sim.h))
colormap('gray')

i=i+ncols; h(i) = subplot(nrows,ncols,i);
Pl.label=[];
Plot_nX(Pl,C);

i=i+ncols; h(i) = subplot(nrows,ncols,i); hold on
kk=1;
for j=1:Nc
    kk=-kk;
    stem(kk*n(:,j),'Marker','none','LineWidth',Pl.sw,'Color',Pl.colors(j,:))
end
axis([Pl.xlims -max(n(:,1)) max(n(:,2))])
ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:1,'YTickLabel',[])
set(gca,'XTick',Pl.XTicks,'XTickLabel',[])
set(gca,'XTickLabel',[])
box off

for j=1:Nc
    
    % plot spatial filter
    i=j+1; h(i) = subplot(nrows,ncols,i);
    imagesc(reshape(Phat{q}.a(:,j),Sim.w,Sim.h)),
    title(I{q}.label)
    if j==1,
        ylab=ylabel([{'Spatial'}; {'Filter'}]);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    else
        set(gca,'XTick',[],'YTick',[])
    end
    
    % plot fluorescence data
    i=i+ncols; h(i) = subplot(nrows,ncols,i);
    Pl.color = 'k';
    Pl.label=[{'Fluorescence'}; {'Projection'}];
    Plot_X(Pl,F*Phat{q}.a(:,j));

    % plot inferred spike trains
    Pl.label = [{'Fast'}; {'Filter'}];
    i=i+ncols; h(i) = subplot(nrows,ncols,i);
    Pl.col(2,:)=[0 0 0];
    Pl.gray=[.5 .5 .5];
    Plot_n_MAP(Pl,I{q}.n(:,j));

    % set xlabel stuff
    subplot(nrows,ncols,i)
    set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
    xlabel('Time (sec)','FontSize',Pl.fs)
    %     linkaxes(h,'x')

    % print fig
    wh=[7 5];   %width and height
    set(fnum,'PaperPosition',[0 11-wh(2) wh]);
    print('-depsc','spatial_multi')
end