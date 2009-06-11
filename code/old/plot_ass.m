clear Pl
nrows   = 3;                                  % set number of rows
ncols   = numel(qs);
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

figure(3), clf
for q=qs
    
    % plot spatial filter
    i=q; h(i) = subplot(nrows,ncols,i);
    imagesc(reshape(Phat{q}.a,Sim.w,Sim.h)),
    if q==1,
        title([{'Optimal Projection'}]);
        Pl.label = 'Fluorescence';
    else
        title([{'SVD Projection'}])
        Pl.label = [];
    end
    ylabel([{'Spatial'}; {'Filter'}])
%     set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
  
    
    % plot fluorescence data
    i=i+ncols; h(i) = subplot(nrows,ncols,i);
    Pl.color = 'k';
    Plot_nX(Pl,F*Phat{q}.a);

    % plot inferred spike trains
    if q==1, Pl.label = [{'Fast'}; {'Filter'}];
    else Pl.label=[]; end
    i=i+ncols; h(i) = subplot(nrows,ncols,i);
    Pl.col(2,:)=[0 0 0];
    Pl.gray=[.5 .5 .5];
    Plot_n_MAP(Pl,I{q}.n);

    % set xlabel stuff
    subplot(nrows,ncols,i)
    set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
    xlabel('Time (sec)','FontSize',Pl.fs)
    %     linkaxes(h,'x')

    % print fig
    wh=[7 5];   %width and height
    set(fnum,'PaperPosition',[0 11-wh(2) wh]);
    print('-depsc','spatial_EM')
end