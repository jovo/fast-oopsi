function [V inf spt volt n_t] =  MakeFigs(V,P)

datasets=V.datasets;
filters=V.filters;
name=V.name;

if ~isempty(datasets)
    dataset = load('~/Research/oopsi/meta-oopsi/data/rafa/adam/2008/Imaging-SNR-Data.mat');
    names = fieldnames(dataset);
end
%%
V.fast_iter_max = 1;
V.save          = 1;
V.plot          = 0;
V.fast_plot     = 1;
V.smc_iter_max  = 5;
V.smc_plot      = 1;
V.N             = 99;
V.est_sig       = 1;
V.est_t         = 0;
P.k_d           = 180;


for i=datasets
    disp(i)
    if i==13

        if ~isfield(V,'T'),     V.T = 400;     end     % # of time steps
        if ~isfield(V,'dt'),    V.dt = 1/30;    end     % # of pixels in ROI
        if ~isfield(V,'fast_iter_max'),    V.fast_iter_max = 1;    end     % # of pixels in ROI

        tau     = 1;                % decay time constant for each cell
        if ~isfield(P,'b'),     P.b     = 0;    end
        if ~isfield(P,'sig'),   P.sig   = 0.4;  end
        if ~isfield(P,'gam'),   P.gam   = 1-V.dt/tau; end
        if ~isfield(P,'lam'),   P.lam   = 2;    end
        if ~isfield(P,'a'),     P.a     = 1;    end

        V.n     = poissrnd(P.lam*V.dt*ones(V.T,1));
        V.n(V.n>1)=1;
        V.C     = filter(1,[1 -P.gam],V.n);         % calcium concentration
        F{i}    = P.a*V.C'+P.b+P.sig*randn(1,V.T);
        spt{i}  = find(V.n);
        n_t{i}  = V.n;
        volt{i} = V.n;

    else
        cc      = dataset.(char(names(i)));
        yy      = cc.Fluorescence;
        V.T     = length(yy);
        f       = detrend(yy);
        nfft    = 2^nextpow2(V.T);
        y       = fft(f,nfft);
        bw      = 10;
        y(1:bw) = 0; y(end-bw+1:end)=0;
        iy      = ifft(y,nfft);
        G       = real(iy(1:V.T));
        ab      = [G'; ones(1,V.T)]'\yy;
        F{i}    = ab(1)*G+ab(2);
        %         F{i}    = z1(cc.Fluorescence);
        %         V.T     = length(F{i});
        %         f       = detrend(F{i});
        %         nfft    = 2^nextpow2(V.T);
        %         y       = fft(f,nfft);
        %         bw      = 10;
        %         y(1:bw) = 0; y(end-bw+1:end)=0;
        %         iy      = ifft(y,nfft);
        %         F{i}    = z1(real(iy(1:V.T)));
        V.dt    = median(diff(cc.FluorescenceTime));
        P.gam   = 1-V.dt/1;
        volt1{i}= cc.chanDev1_ai0_VoltageCh1;
        volt{i} = interp1(cc.time,volt1{i},cc.FluorescenceTime);
        temp    = GetSpikeTimes(volt1{i},0.7);
        spt{i}=[];
        for t=1:length(temp)
            spt{i} = [spt{i} find(cc.time(temp(t))<cc.FluorescenceTime,1)];
        end
        n_t{i}  = zeros(size(F{i}));
        n_t{i}(spt{i}) = 1;
        V.n = n_t{i};
        V.name = [name num2str(i)];
    end
    for j=filters
        switch j
            case 0 % fast init smc
                V.fast_iter_max = 2;
                V.smc_iter_max  = 3;
                %                 V.T=min(length(F{i}),1000);
                [fast smc] = run_oopsi(F{i}(1:V.T),V,P);
                V.fast_iter_max = 1;
                inf{i}.fast = fast;
                inf{i}.smc = smc;
                PP = fast.P;
            case 1 % fast with params
                V.fast_poiss=0;
                V.fast_nonlin=0;
                V.test=0;
                [inf{i}.fast.n inf{i}.fast.P] = fast_oopsi(F{i},V,P);
            case 1.1
                V.fast_poiss=0;
                V.fast_nonlin=0;
                V.test=1;
                [inf{i}.fast1 PP] = fast_oopsi(F{i},V,P);
            case 1.5 % fast est params
                V.fast_iter_max=1;
                V.fast_plot=1;
                V.fast_poiss=0;
                V.fast_nonlin=0;
                V.b_est=0;
                V.sig_est=1;
                [inf{i}.fast.n inf{i}.fast.P] = fast_oopsi(F{i},V);
            case 2 % nonlin
                V.fast_poiss=0;
                V.fast_nonlin=1;
                V.gauss_n = inf{i}.fast + 0.0001;
                C       = filter(1,[1 -P.gam],V.gauss_n)';         % calcium concentration
                C1      = [C./(C+P.k_d); ones(1,V.T)];         % for brevity
                ab      = C1'\F{i}';                           % estimate scalse and offset
                PN      = PP;
                PN.a     = ab(1);                               % fluorescence scale
                PN.b     = ab(2);                               % fluorescence offset
                inf{i}.nonlin = fast_oopsi(F{i},V,PN);
            case 3 % poisson
                V.fast_poiss=1;
                V.fast_nonlin=0;
                inf{i}.poisson = fast_oopsi(F{i},V,P);
            case 4 % Wiener known params
                PP.lam = P.sig;
                %                 PP.lam = sum(inf{i}.fast/max(inf{i}.fast))/(V.T*V.dt);
                PP=inf{i}.fast.P;
                %                 PP.sig = inf{i}.fast.P.sig;
                %                 PP.gam = inf{i}.fast.P.gam
                [inf{i}.Wiener.n inf{i}.Wiener.P] = wiener_oopsi(F{i},V.dt,PP);
            case 4.5 % Wiener unknown params
                PP = inf{i}.fast.P;
                PP.sig = mad(F{i},1)*1.4826;
                PP.lam = PP.sig;
                [inf{i}.Wiener.n inf{i}.Wiener.P] = wiener_oopsi(F{i},V.dt,PP);
            case 5 % smc
                [M P V] = smc_oopsi(F{i},V,P);
                inf{i}.smc = M.nbar;
        end
    end
    V.name = [name num2str(i)];
    if V.save==1, save(['../../data/' V.name],'F','inf','spt','volt','n_t','V'); end
end


%%
% datasets=12;
% V.name=['fast_smc_vitro' num2str(datasets)];
% load(['../../data/' V.name])

for j=datasets
    V.name_fig = ['../../figs/' V.name];                                 % filename for figure
    fig     = figure(j); clf,
    V.T     = length(F{j});
    names   = fieldnames(inf{j});
    nrows   = 2 + length(names);
    gray    = 0.4*[1 1 1];              % define gray color
    inter   = 'tex';                    % interpreter for axis labels
    fs      = 14;                       % font size
    ms      = 5;                        % marker size for real spike
    sw      = 2;                        % spike width
    lw      = 1;                        % line width
    if ~isfield(V,'xlims') xlims   = [200 V.T]; else xlims = V.xlims; end
    xticks  = xlims(1):1/V.dt:xlims(2);             % XTick positions
    skip    = round(length(xticks)/5);
    xticks  = xticks(1:skip:end);
    tvec_o  = xlims(1):xlims(2);        % only plot stuff within xlims
    maxn    = max(V.n);

    % plot fluorescence data
    i=1; h(i)=subplot(nrows,1,i); hold on
    Ftemp=z1(F{j}(tvec_o))*maxn;
    Ftemp=Ftemp-mean(Ftemp);
    maxF=max(Ftemp);
    minF=min(Ftemp);
    miny=min(minF,floor(minF)/2);
    maxy=min(maxF,ceil(maxF)/2);
    hold on
    stem(spt{j},n_t{j}(spt{j})*(maxF),'Marker','+','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
    plot(zeros(V.T,1),'w')
    plot(tvec_o,Ftemp,'-k','LineWidth',lw);
    ylab=ylabel([{'fluorescence'}; {'(%dF/F)    '}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    XTicks=[floor(minF):.5:ceil(maxF)];
    set(gca,'YTick',XTicks,'YTickLabel',XTicks*100)
    set(gca,'XTick',xticks,'XTickLabel',[])
    axis([xlims miny max(maxy,maxF)])

    % plot voltage data
    i=i+1; h(i)=subplot(nrows,1,i); hold on
    volt{j}(spt{j})=max(volt{j});
    if datasets==13
        stem(V.n,'Marker','none','LineWidth',sw,'Color','k')
        ylab=ylabel([{'spike train'}; {'(# spikes)'}],'Interpreter',inter,'FontSize',fs);
        set(gca,'YTick',[0:maxn],'YTickLabel',[0:maxn])
    else
        plot(z1(volt{j}),'-k','LineWidth',lw);
        ylab=ylabel([{'voltage'}; {'(mV)  '}],'Interpreter',inter,'FontSize',fs);
        set(gca,'YTick',[0:maxn],'YTickLabel',round([min(volt{j}) max(volt{j})]))
    end
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'XTick',xticks,'XTickLabel',[])
    axis([xlims 0 maxn])

    for k=1:length(names)
        i=i+1; h(i)=subplot(nrows,1,i); hold on,
        n = inf{j}.(char(names(k))).n;
        n = n/max(abs(n(tvec_o)))*maxn;

        nneg=find(n<0);
        stem(nneg,n(nneg),'Marker','none','LineWidth',sw,'Color',gray)

        npos=find(n>0);
        stem(npos,n(npos),'Marker','none','LineWidth',sw,'Color','k')

        stem(spt{j},n_t{j}(spt{j}),'Marker','+','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
        axis([xlims floor(min(n(tvec_o))) max(n(tvec_o))])
        hold off,
        if strcmp(char(names(k)),'fast')
            ylab=ylabel([[char(names(k)) ' filter']; {'(a.u.)   '}],'Interpreter',inter,'FontSize',fs);
        else
            ylab=ylabel([[char(names(k)) ' filter']; {'(a.u.)     '}],'Interpreter',inter,'FontSize',fs);
        end
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
        set(gca,'YTick',[floor(min(n)):ceil(maxn)],'YTickLabel',[floor(min(n)):ceil(maxn)])
        set(gca,'XTick',xticks,'XTickLabel',[])
        box off
    end

    % label last subplot
    set(gca,'XTick',xticks,'XTickLabel',(xticks-xticks(1))*V.dt)
    xlabel('time (sec)','FontSize',fs)
    linkaxes(h,'x')

    % print fig
    if V.save==1
        wh=[7 5];   %width and height
        set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
        print('-depsc',V.name_fig)
        print('-dpdf',V.name_fig)
        saveas(fig,V.name_fig)
    end
end