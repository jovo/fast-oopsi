% SI_ANALYSIS Test Josh Vogelstein's spike inference algorithms
% given a standard data structure. Plot electrophysiological data
% and calcium transients against with the inferred spike times from
% both the fast and slow spike inference techniques (oopsi).
%
% tamachado

clear; clear global; clc

dataset = load('~/Research/oopsi/meta-oopsi/data/rafa/adam/2008/Imaging-SNR-Data.mat');
% load '~/Research/oopsi/meta-oopsi/docs/posters/SFN09/oopsi.mat';
names = fieldnames(dataset);
%%
V.wiener_do     = 1;
V.fast_do       = 1;
V.smc_do        = 0;
V.save          = 0;
V.plot          = 0;
V.fast_plot     = 1;
V.fast_iter_max = 1;
V.smc_iter_max  = 25;
V.N             = 99;
V.sig_est       = 1;
P.lam           = 10;
P.k_d           = 180;
datasets        = 1:12;

for i=datasets
    cc      = dataset.(char(names(i)));
    F{i}    = z1(cc.Fluorescence);
    V.dt    = median(diff(cc.FluorescenceTime));
    P.gam   = 1-V.dt/1;
    filters = [1 2 3 4];
    for j=filters
        switch j
            case 1
                V.fast_poiss=0;
                V.fast_nonlin=0;
                [inf{i}.vanilla P] = fast_oopsi(F{i},V,P);
            case 2
                V.fast_poiss=1;
                V.fast_nonlin=0;
                inf{i}.poisson = fast_oopsi(F{i},V,P);
            case 3
                V.fast_poiss=0;
                V.fast_nonlin=1;
                V.gauss_n = inf{i}.vanilla + 0.0001;
                inf{i}.nonlin = fast_oopsi(F{i},V,P);
            case 4
                inf{i}.wiener = wiener_oopsi(F{i},V.dt,P);
            case 5
                [M P V] = smc_oopsi(F{i},V,P);
                inf{i}.smc = M.nbar;
        end
    end
    volt1{i}= cc.chanDev1_ai0_VoltageCh1;
    volt{i} = interp1(cc.time,volt1{i},cc.FluorescenceTime);
    temp    = GetSpikeTimes(volt1{i},0.7);
    spt{i}=[];
    for t=1:length(temp)
        spt{i} = [spt{i} find(cc.time(temp(t))<cc.FluorescenceTime,1)];
    end
    n_t{i}  = zeros(size(inf{i}.vanilla));
    n_t{i}(spt{i}) = 1;
    save(['../../data/adam'],'F','inf','spt','volt')
end


%%
for j=datasets
    V.name_fig = ['../../figs/adam' num2str(j)];                                 % filename for figure
    fig     = figure(j); clf,
    V.T     = length(F{j});
    names   = fieldnames(inf{j});
    nrows   = 2 + 2; %length(names);
    gray    = [.75 .75 .75];            % define gray color
    inter   = 'tex';                    % interpreter for axis labels
    xlims   = [5 V.T-15];                  % xmin and xmax for current plot
    fs      = 14;                       % font size
    ms      = 5;                        % marker size for real spike
    sw      = 2;                        % spike width
    lw      = 2;                        % line width
    xticks  = xlims(1):1/V.dt:xlims(2);             % XTick positions
    skip    = round(length(xticks)/5);
    xticks  = xticks(1:skip:end);
    tvec_o  = xlims(1):xlims(2);        % only plot stuff within xlims

    % plot fluorescence data
    i=1; h(i)=subplot(nrows,1,i); hold on
    plot(tvec_o,z1(F{j}(tvec_o)),'-k','LineWidth',lw);
    ylab=ylabel([{'fluorescence'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[],'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    axis([xlims 0 1.1])

    % plot voltage data
    i=i+1; h(i)=subplot(nrows,1,i); hold on
    volt{j}(spt{j})=max(volt{j});
    plot(z1(volt{j}),'-k','LineWidth',lw);
    ylab=ylabel([{'voltage'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[],'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    axis([xlims 0 1.1])

    for k=filters
        i=i+1; h(i)=subplot(nrows,1,i); hold on,
        n = inf{j}.(char(names(k)));
        n = n/max(abs(n(tvec_o)));
        
        nneg=find(n<0);
        stem(nneg,n(nneg),'Marker','none','LineWidth',sw,'Color',gray)

        npos=find(n>0);
        stem(npos,n(npos),'Marker','none','LineWidth',sw,'Color','k')

        stem(spt{j},n_t{j}(spt{j}),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
        axis([xlims min(n) 1.1])
        hold off,
        ylab=ylabel([{char(names(k))}; {'filter'}],'Interpreter',inter,'FontSize',fs);
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
        set(gca,'YTick',0:2,'YTickLabel',[])
        set(gca,'XTick',xticks,'XTickLabel',[])
        box off
    end

    % label last subplot
    set(gca,'XTick',xticks,'XTickLabel',(xticks-xticks(1))*V.dt)
    xlabel('Time (sec)','FontSize',fs)
    linkaxes(h,'x')

    % print fig
    wh=[7 5];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    print('-depsc',V.name_fig)
    print('-dpdf',V.name_fig)
    saveas(fig,V.name_fig)
end