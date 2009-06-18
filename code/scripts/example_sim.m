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
Ntrials=1;
ks=2.^(0:4);
for kk=2%1:length(ks)
    disp(kk)
    for tt=1:Ntrials
        disp(tt)
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
        Sim.T       = 500;                              % # of time steps
        Sim.dt      = 0.005;                            % time step size
        Sim.MaxIter = 0;                                % # iterations of EM to estimate params
        Sim.Np      = Npixs;                            % # of pixels in each image
        Sim.w       = width;                            % width of frame (pixels)
        Sim.h       = height;                           % height of frame (pixels)
        Sim.Nc      = Nc;                               % # cells
        Sim.plot    = 0;                                % whether to plot filter with each iteration

        lam         = [10; Sim.T];
        sigs        = [1/8 4]*ks(kk);
        moda        = (sin(linspace(0,10*pi,Sim.T-1))+1)/2;
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
            tau     = [.1 .5]; %round(100*rand(Sim.Nc,1))/100+0.05;   % decay time constant for each cell
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


        end


        %% end) plot results
        if tt==Ntrials && kk==2
            clear Pl
            nrows   = 2+Nc;                                  % set number of rows
            ncols   = 2;
            h       = zeros(nrows,1);
            Pl.xlims= [5 Sim.T-101];                            % time steps to plot
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

                % plot fast spike trains
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
                    plot(ass(Sim.T/2+(Pl.xlims(1):Pl.xlims(2))),'Color',Pl.gray,'LineWidth',Pl.lw);
                    ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs);
                    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
                    set(gca,'YTick',[],'YTickLabel',[])
                    set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
                    X=[I{1,q}.n D{tt,q}.n];
                    axis([Pl.xlims-Pl.xlims(1) min(X(:)) max(ass(:))])
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

                % plot wiener spike trains
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
                    plot(ass(Sim.T/2+(Pl.xlims(1):Pl.xlims(2))),'Color',Pl.gray,'LineWidth',Pl.lw);
                    ylab=ylabel(Pl.label,'Interpreter',Pl.inter,'FontSize',Pl.fs);
                    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
                    set(gca,'YTick',[],'YTickLabel',[])
                    set(gca,'XTick',Pl.XTicks,'XTickLabel',[],'FontSize',Pl.fs)
                    X=[I{2,q}.n D{tt,q}.n];
                    axis([Pl.xlims-Pl.xlims(1) min(X(:)) max(ass(:))])
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

            % print fig
            wh=[7 5];   %width and height
            DirName = '../../figs/';
            FileName = 'example_sim';
            PrintFig(wh,DirName,FileName);
        end
    end


end


