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

clear,
clc
Ntrials=5;
ks=2.^(0:4);
sig0=[1/16 1/2];
for kk=1:length(ks)
    disp(['noise ' num2str(kk)])
    for tt=1:Ntrials
        disp(['trial ' num2str(tt)])
        % 1) generate spatial filters

        % % stuff required for each spatial filter
        Nc      = 1;                                % # of cells in the ROI
        neur_w  = 1;                               % width per neuron
        width   = 1;                               % width of frame (pixels)
        height  = Nc*neur_w;                        % height of frame (pixels)
        Npixs   = width*height;                     % # pixels in ROI
        % x1      = linspace(-5,5,height);
        % x2      = linspace(-5,5,width);
        % [X1,X2] = meshgrid(x1,x2);
        % g1      = zeros(Npixs,Nc);
        % g2      = 0*g1;
        % Sigma1  = diag([1,1])*1;                    % var of positive gaussian
        % Sigma2  = diag([1,1])*2;                    % var of negative gaussian
        % mu      = [1 1]'*linspace(-2,2,Nc);         % means of gaussians for each cell (distributed across pixel space)
        % w       = Nc:-1:1;                          % weights of each filter
        %
        % % spatial filter
        % for i=1:Nc
        %     g1(:,i)  = w(i)*mvnpdf([X1(:) X2(:)],mu(:,i)',Sigma1);
        %     g2(:,i)  = w(i)*mvnpdf([X1(:) X2(:)],mu(:,i)',Sigma2);
        % end
        % a_b = sum(g1-g2,2);

        % 2) set simulation metadata
        Sim.T       = 5000;                              % # of time steps
        Sim.dt      = 0.05;                            % time step size
        Sim.MaxIter = 0;                                % # iterations of EM to estimate params
        Sim.Np      = Npixs;                            % # of pixels in each image
        Sim.w       = width;                            % width of frame (pixels)
        Sim.h       = height;                           % height of frame (pixels)
        Sim.Nc      = Nc;                               % # cells
        Sim.plot    = 0;                                % whether to plot filter with each iteration

        lam         = [1; 10];
        sigs        = sig0*ks(kk);
        moda        = (sin(linspace(0,10*pi,Sim.T-1))+1);
        %         moda        = mean(moda)+0*moda;
        qs          = 1:2;

        for q=qs
            % 3) initialize params
            P.a     = 1;
            % for i=1:Sim.Nc
            %     P.a(:,i)=g1(:,i)-g2(:,i);
            % end
            P.b     = 0;                           % baseline is zero
            % P.b     = 0*P.a(:,1)+1;                           % baseline is zero

            P.sig   = sigs(q);                                 % stan dev of noise (indep for each pixel)
            C_0     = 0;                                    % initial calcium
            tau     = 1; %round(100*rand(Sim.Nc,1))/100+0.05;   % decay time constant for each cell
            P.gam   = 1-Sim.dt./tau(1:Sim.Nc);
            P.lam   = lam(q);%round(10*rand(Sim.Nc,1))+5;           % rate-ish, ie, lam*dt=# spikes per second

            % 3) simulate data
            n=zeros(Sim.T,Sim.Nc);
            C=n;
            for i=1:Sim.Nc
                n(1,i)      = C_0;
                n(2:end,i)  = poissrnd(P.lam(i)*Sim.dt*moda);    % simulate spike train
                if q==1, n(n>1)=1; end
                C(:,i)      = filter(1,[1 -P.gam(i)],n(:,i));               % calcium concentration
            end
            Z = 0*n(:,1);
            F = C*P.a' + (1+Z)*P.b'+P.sig*randn(Sim.T,Npixs);               % fluorescence

            D{tt,q}.n=n; D{tt,q}.C=C; D{tt,q}.F=F;

            %% infer spikes
            GG=D{tt,q}.F; Tim=Sim;
            Phat{tt,q}=P;
            I{tt,q}.label='True Filter';
            display(I{tt,q}.label)
            tic
            I{1,q}.n = FOOPSI2_59(GG,Phat{tt,q},Tim);
            I{1,q}.time = toc;

            tic
            I{2,q}.n = WienerFilt1_2(F,Sim.dt,P);
            I{2,q}.time = toc;

            I{3,q}.n = I{2,q}.n; I{3}.n(I{3,q}.n<0)=0;
            I{3,q}.time = toc;

            tic
            I{4,q}.n = [diff(F); 0];
            I{4,q}.time = toc;

            a = 75;
            gaus = exp(-(-Sim.T/2:Sim.T/2).^2);

            tic
            I{5,q}.n = conv2(exp(-(linspace(-a,a,Sim.T)).^2)',I{1,q}.n,'same');
            I{5,q}.time = toc;

            tic
            I{6,q}.n = conv2(exp(-(linspace(-a,a,Sim.T)).^2)',I{2,q}.n,'same');
            I{6,q}.time = toc;

            %% compute statistics
            for i=1:6

                if q==1,
                    I{i,q}.n=I{i,q}.n/max(abs(I{i,q}.n));
                else
                    %                 I{i,q}.n=I{i,q}.n/(max(abs(I{i,q}.n))*max(D{tt,q}.n));
                end

                xc=corrcoef(I{i,q}.n,D{tt,q}.n);
                cc{q}(tt,i)=xc(2);

                snr{q}(tt,i)=mean(I{i,q}.n(D{tt,q}.n==1).^2)/mean(I{i,q}.n(D{tt,q}.n==0).^2);

                roc{q,tt,i} = roc2(I{i,q}.n, D{tt,q}.n);

                mse{q}(tt,i)=mean((I{i,q}.n-D{tt,q}.n).^2);

                time{q}(tt,i)=I{i,q}.time;
            end

        end

        %% end) plot results
        if tt==11 && kk==2
            clear Pl
            nrows   = 2+Nc;                                  % set number of rows
            ncols   = 2;
            h       = zeros(nrows,1);
            Pl.xlims= [3001 3300];                            % time steps to plot
            Pl.nticks=5;                                    % number of ticks along x-axis
            Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
            Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
            Pl.vs   = 4;
            Pl.colors(1,:) = [0 0 0];
            Pl.colors(2,:) = Pl.gray;
            Pl.colors(3,:) = [.5 0 0];
            Pl.Nc   = Sim.Nc;
            fnum = figure(1); clf,
            Pl.interp = 'latex';

            for q=qs
                Pl.n    = double(D{tt,q}.n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting

                % plot fluorescence data
                i=q; h(i) = subplot(nrows,ncols,i);
                if q==1,
                    Pl.label = [{'$\mathbf{F}$'}];
                    title('Slow Firing Rate')
                else
                    Pl.label = [];
                    title('Fast Firing Rate')
                end
                Pl.color = 'k';
                Plot_nX(Pl,D{tt,q}.F*Phat{tt,q}.a);
                %     title(I{q}.label)
                axis('tight')

                % plot FAND filter
                i=i+2; h(i) = subplot(nrows,ncols,i);
                if q==1,
                    Pl.label = [{'$\mathbf{n}_{FAND}$'}];
                    Pl.j=1;
                    Plot_n_MAPs(Pl,I{1,q}.n);
                else
                    Pl.label = [];
                    hold on
                    a=75;
                    gaus = exp(-(-Sim.T/2:Sim.T/2).^2);
                    ass=conv(exp(-(linspace(-a,a,Sim.T)).^2),I{1,q}.n);
                    plot(ass(Sim.T/2+(Pl.xlims(1):Pl.xlims(2))),'Color','k','LineWidth',Pl.lw);
                    ass=conv(exp(-(linspace(-a,a,Sim.T)).^2),D{tt,q}.n);
                    plot(ass(Sim.T/2+(Pl.xlims(1):Pl.xlims(2))),'Color',Pl.gray,'LineWidth',1);
                    ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs);
                    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
                    set(gca,'YTick',[],'YTickLabel',[])
                    set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
                    X=[I{1,q}.n D{tt,q}.n];
                    axis([Pl.xlims-Pl.xlims(1) min(X(:)) max(ass(:))])
                    axis('tight')
                    box off
                end

                % set ylabel stuff
                if q==1
                    Pl.YTicks=0:max(D{q}.n):max(D{q}.n);
                else
                    Pl.YTicks=0:round(max(ass)/2):round(max(ass));
                end
                Pl.YTickLabels=Pl.YTicks;
                set(gca,'YTick',Pl.YTicks,'YTickLabel',Pl.YTicks, 'FontSize',10)

                % plot wiener filter trains
                i=i+2; h(i) = subplot(nrows,ncols,i);
                if q==1,
                    Pl.label = [{'$\mathbf{n}_{Wiener}$'}];
                    Pl.j=1;
                    Plot_n_MAPs(Pl,I{q+2}.n);
                else
                    Pl.label = [];
                    hold on
                    a=75;
                    gaus = exp(-(-Sim.T/2:Sim.T/2).^2);
                    ass=conv(exp(-(linspace(-a,a,Sim.T)).^2),I{2,q}.n);
                    plot(ass(Sim.T/2+(Pl.xlims(1):Pl.xlims(2))),'Color','k','LineWidth',Pl.lw);
                    ass=conv(exp(-(linspace(-a,a,Sim.T)).^2),D{tt,q}.n);
                    plot(ass(Sim.T/2+(Pl.xlims(1):Pl.xlims(2))),'Color',Pl.gray,'LineWidth',1);
                    ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs);
                    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
                    set(gca,'YTick',[],'YTickLabel',[])
                    set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
                    X=[I{2,q}.n D{tt,q}.n];
                    axis([Pl.xlims-Pl.xlims(1) min(X(:)) max(ass(:))])
                    axis('tight')
                    box off
                end

                % set ylabel stuff
                if q==1
                    Pl.YTicks=0:max(D{q}.n):max(D{q}.n);
                else
                    Pl.YTicks=0:round(max(ass)/2):round(max(ass));
                end
                Pl.YTickLabels=Pl.YTicks;
                set(gca,'YTick',Pl.YTicks,'YTickLabel',Pl.YTicks, 'FontSize',10)

                subplot(nrows,ncols,i)
                set(gca,'XTick',Pl.XTicks,'XTickLabel',Pl.XTicks*Sim.dt,'FontSize',Pl.fs)
                xlabel('Time (sec)','FontSize',Pl.fs)
                %     linkaxes(h,'x')

            end

            %             % print fig
            %             wh=[7 5];   %width and height
            %             DirName = '../../figs/';
            %             FileName = 'wiener';
            %             PrintFig(wh,DirName,FileName);
        end
    end


    %% make errorbar stats

    for tt=1:Ntrials
        for i=1:6
            auc.tn(tt,i)=sum(roc{1,tt,i}.tn);
            auc.tp(tt,i)=sum(roc{1,tt,i}.tp);
        end
    end

    for q=qs
        stats.cc{q}.means(kk,:) = mean(cc{q});
        stats.cc{q}.std(kk,:)   = std(cc{q});

        stats.mse{q}.means(kk,:) = mean(mse{q});
        stats.mse{q}.std(kk,:)   = std(mse{q});

        stats.snr{q}.means(kk,:) = mean(snr{q});
        stats.snr{q}.std(kk,:)   = std(snr{q});

        stats.auc{q}.means(kk,:) = mean(auc.tp);
        stats.auc{q}.std(kk,:)   = std(auc.tp);

        stats.time{q}.means(kk,:) = mean(time{q});
        stats.time{q}.std(kk,:)   = std(time{q});

    end

    %% make table

    %     if kk==2
    %         for q=qs
    %             y=zeros(4,10);
    %
    %             y(:,1)=mean(cc{q});
    %             y(:,2)=std(cc{q});
    %
    %             y(:,3)=mean(mse{q});
    %             y(:,4)=std(mse{q});
    %
    %             y(:,5)=mean(snr{q});
    %             y(:,6)=std(snr{q});
    %
    %             y(:,7)=mean(auc.tp);
    %             y(:,8)=std(auc.tp);
    %
    %             y(:,9)=mean(time{q});
    %             y(:,10)=std(time{q});
    %
    %             if q==1
    %                 fid = fopen([DirName 'tab_slow.tex'],'w');
    %             elseif q==2
    %                 fid = fopen([DirName 'tab_fast.tex'],'w');
    %                 y(:,[5:8])=NaN;
    %             end
    %             fprintf(fid,'Filter & $R_{\\hbn}$ & $MSE$ & $SNR$ & $AUC$ & Time \\\\ \\hline \n');
    %
    %             for i=1:4
    %                 if i==1
    %                     fprintf(fid,'FAND');
    %                 elseif i==2
    %                     fprintf(fid,'Wiener');
    %                 elseif i==3
    %                     fprintf(fid,'$[$Wiener$]_+$');
    %                 elseif i==4
    %                     fprintf(fid,'dF/F');
    %                 end
    %                 fprintf(fid,'& %2.3f (%2.3f) & %2.3f (%2.3f) & %2.3f (%2.3f) & %2.3f (%2.3f) & %2.3f (%2.3f)',y(i,:));
    %
    %                 if i<4, fprintf(fid,' \\\\ \n'); end
    %             end
    %             status = fclose(fid);
    %         end
    %     end

end
save('../../data/stats_fig.mat')

%% make fig
load('../../data/stats_fig.mat')

Pl.fs = 18;

figure(2), clf
nrows=2; ncols=2;
for q=qs

    if q==1
        subplot(nrows,ncols,1),
        errorbar(stats.snr{q}.means(:,1:4),stats.snr{q}.std(:,1:4),'LineWidth',2)
        ymin=min(min(stats.snr{q}.means(:,1:4)+stats.snr{q}.std(:,1:4)));
        ymax=max(max(stats.snr{q}.means(:,1:4)+stats.snr{q}.std(:,1:4)));
        axis([1 5  ymin*.9 ymax])
        set(gca,'YScale','log','YTick',10.^(-5:10),'XTickLabel',[],'FontSize',Pl.fs);
        set(gca,'YTick',10.^(0:1:5),'YTickLabel',10.^(0:1:5))
        ylab=ylabel('eSNR','Interpreter','none');
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',Pl.fs)
        title('Sparse Spiking','FontSize',Pl.fs)
    end

    %     % plot corrcoef
    %     subplot(nrows,ncols,1+(q-1)*ncols),
    %     errorbar(stats.cc{q}.means(:,1:4),stats.cc{q}.std(:,1:4))
    %     if q==1,
    %         title('$\rho$','Interpreter','latex'),
    %         ylab=ylabel([{'Sparsely'}; {'Spiking'}; {'Neuron'}]);
    %     elseif q==2
    %         ylab=ylabel([{'Fast'}; {'Spiking'}; {'Neuron'}]);
    %     end
    %     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',12)
    %
    %     axis([1 5 0 1])
    %     set(gca,'YScale','linear','YTick',[0 .5 1],'XTickLabel',[]);

    % plot mse
    if q==2
        subplot(nrows,ncols,ncols),
        errorbar(stats.mse{q}.means(:,1:4),stats.mse{q}.std(:,1:4),'LineWidth',2)
        if q==1, title('MSE','Interpreter','none'), end
        ymin=min(min(stats.mse{q}.means(:,1:4)+stats.mse{q}.std(:,1:4)));
        ymax=max(max(stats.mse{q}.means(:,1:4)+stats.mse{q}.std(:,1:4)));
        axis([1 5 ymin*.5 ymax])
        set(gca,'YScale','log','YTick',10.^(-5:10),'XTickLabel',[]);
        set(gca,'YTick',10.^(0:1:5),'YTickLabel',10.^(0:1:5))
        if q==2,
            set(gca,'XTick',[1:5], 'XTickLabel',[{'4'}; {''}; {'16'}; {''}; {'32'}],'FontSize',Pl.fs);
            xlabel('$\sigma$','Interpreter','latex');
            ylab=ylabel('MSE','Interpreter','none');
            set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontSize',Pl.fs)
        end
        title('Fast Spiking','FontSize',Pl.fs)
    end
    %     if q==1
    %         subplot(nrows,ncols,4+(q-1)*ncols),
    %         errorbar(stats.auc{q}.means,stats.auc{q}.std)
    %         ymin=min(stats.auc{q}.means(:)+stats.auc{q}.std(:));
    %         ymax=max(stats.auc{q}.means(:)+stats.auc{q}.std(:));
    %         axis([1 5  0 10])
    %         set(gca,'YScale','linear','YTick',[0 5 10],'YTickLabel',[0 5 10])
    %         title('AUC','Interpreter','none')
    %         xlabel('$\sigma$','Interpreter','latex')
    %     end

    %     subplot(nrows,ncols,5+(q-1)*ncols),
    %     errorbar(stats.time{q}.means,stats.time{q}.std)
    %     ymin=min(stats.time{q}.means(:)+stats.time{q}.std(:));
    %     ymax=max(stats.time{q}.means(:)+stats.time{q}.std(:));
    %     axis([1 5 ymin*.5 ymax])
    %     set(gca,'YScale','log','YTick',10.^(-5:10),'XTickLabel',[]);
    %     if q==1, title('Time (sec)','Interpreter','none'), end
    if q==1
        set(gca,'XTick',[1:2:5], 'XTickLabel',[{'1/16'}; {'1/4'}; {'1'}]);
    else
        set(gca,'XTick',[1:2:5], 'XTickLabel',sig0(q)*ks(1:2:5));
    end
    xlabel('$\sigma$','Interpreter','latex')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot snr subfig
load('../../data/snr.mat')
nrows=2; ncols=2;
clear Pl
clc
Pl.xlims= [5 Sim.T-101];                            % time steps to plot
Pl.nticks=5;                                    % number of ticks along x-axis
Pl.n    = double(n); Pl.n(Pl.n==0)=NaN;         % store spike train for plotting
Pl      = PlotParams(Pl);                       % generate a number of other parameters for plotting
Pl.vs   = 4;
Pl.colors(1,:) = [0 0 0];
Pl.colors(2,:) = Pl.gray;
Pl.colors(3,:) = [.5 0 0];
Pl.Nc   = Sim.Nc;
Pl.interp = 'latex';
Pl.fs = 18;

subplot(nrows,ncols,ncols+1)
h=errorbar(mean_snr,std_snr);
set(gca,'YScale','log') %,'YTick',10.^(-5:10),'XTickLabel',[]);
set(h,'LineWidth',2)
ymax=max(mean_snr(:)+std_snr(:));
ymin=max(.8,min(mean_snr(:)-std_snr(:)));
axis([.9 5.1 ymin ymax*1.1])
% set(gca,'YTick',logspace(log(ymin),log(ymax),5),'YTickLabel',0:5)
set(gca,'XTick',1:2:5,'XTickLabel',lam(1:2:end),'FontSize',Pl.fs)
set(gca,'YTick',10.^(0:1:5),'YTickLabel',10.^(0:1:5))
xlabel('$\lambda$','FontSize',Pl.fs,'Interpreter','latex')
ylab=ylabel([{'eSNR'}],'Interpreter','none','FontSize',Pl.fs,'Color','k');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot time subfig
load('../../data/time_stuff.mat')
nrows=2; ncols=2; Pl.fs=18;
% annotation('line',[0 1],.47*ones(2,1))
% annotation('line',.5*ones(2,1),[0 .5])


subplot(nrows,ncols,nrows*ncols)
% h=errorbar(repmat([1:5],4,1)',mean_snr',mean_snr'/100,std_snr');

errorbar(repmat([1:8]',1,4),mean_time(1:4,:)',mean_time(1:4,:)'/100, std_time(1:4,:)','LineWidth',2)
set(gca,'YScale','log') %,'YTick',10.^(-5:10),'XTickLabel',[]);

ymax=max(mean_time(:)+std_time(:));
ymin=max(10^-5,min(mean_time(:)-std_time(:)));
axis([.9 8.1 ymin ymax])
% set(gca,'YTick',linspace(mink,maxk,5),'YTickLabel',0:5)
set(gca,'XTick',[1 4 8],'XTickLabel',200*2.^([1 4 8]))
yticks=10.^(-4:2:5);
set(gca,'YTick',yticks,'YTickLabel',[{'1e-4'}; {'0.01'}; {'1'}],'FontSize',Pl.fs)
xlabel('# of frames','FontSize',Pl.fs)
ylab=ylabel(['Time (sec)'],'Interpreter','none','FontSize',Pl.fs);
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')


% print fig
wh=[7 5];   %width and height
DirName = '../../figs/';
FileName = 'stats';
PrintFig(wh,DirName,FileName);
