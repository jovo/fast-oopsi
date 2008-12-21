function Enew = GOOPSI_M_v1_0(Sim,S,M,E,F)
% this function finds the mle of the parameters
%
% Input---
% Sim:  simulation parameters
% R:    real data
% S:    simulation results
% M:    moments and sufficient stats
% E:    old parameter estimates
%
% Output is 'Enew', a structure with a field for each parameter, plus some
% additional fields for various likelihoods

Enew = E;   % initialize parameters
lik = [];   % initialize likelihood

if Sim.n_params == true
    % MLE for spike rate parameters: baseline (b), linear filter (k), and spike history weights (omega)
    fprintf('\nestimating spike rate params\n')
    RateParams=E.k;
    sp      = S.n==1;                                               % find (particles,time step) pairs that spike
    nosp    = S.n==0;                                               % don't spike
    x       = repmat(Sim.x,1,Sim.N);                                % generate matrix for gradinent
    zeroy   = zeros(Sim.N,Sim.T);                                   % make matrix of zeros for evaluating lik

    if Sim.h_params == true
        if Sim.M>0                                                  % if spike history terms are present
            RateParams=[RateParams; E.omega];                       % also estimate omega
            for i=1:Sim.M                                           % and modify stimulus matrix for gradient
                x(Sim.StimDim+i,:)=reshape(S.h(:,:,i),1,Sim.N*Sim.T);
            end
        end

        options         = optimset('Display','off','GradObj','off');% use gradient
        [bko lik_r]  = fminunc(@f_bko,RateParams,options);     % find MLE
        Enew.k          = bko(1:end-Sim.M);                         % set new parameter estimes
        if Sim.M>0                                                  % for omega too
            Enew.omega = bko(end-Sim.M+1:end);
        end

    else
        if Sim.M>0                                                  % if spike history terms are present
            for i=1:Sim.M                                           % and modify stimulus matrix for gradient
                x(Sim.StimDim+i,:)=reshape(S.h(:,:,i),1,Sim.N*Sim.T);
            end
        end

        options         = optimset('Display','off','GradObj','off');% use gradient
        [bk lik_r]   = fminunc(@f_bk,RateParams,options);      % find MLE
        Enew.k          = bk(1:end);                                % set new parameter estimes
    end
    Enew.lik_r   = -lik_r;
    lik = [lik Enew.lik_r];
end

    function [lik dlik]= f_bko(RateParams)                          % get lik and grad

        xk      = RateParams(1:end-Sim.M)'*Sim.x;                   % filtered stimulus
        hs      = zeroy;                                            % incorporate spike history terms
        for l=1:Sim.M
            hs  = hs+RateParams(end-Sim.M+l)*S.h(:,:,l);
        end
        s       = repmat(xk,Sim.N,1) + hs;

        f_kdt   = exp(s)*Sim.dt;                                    % shorthand
        lik     = -sum(S.w_b(sp).*log(1-exp(-f_kdt(sp))))...        % liklihood
            +sum(S.w_b(nosp).*f_kdt(nosp));

        if nargout > 1                                              % if gradobj=on
            ef      = exp(f_kdt);                                   % shorthand
            dlik      = x(:,sp)*(S.w_b(sp).*f_kdt(sp)./( ef(sp)-1))... %gradient of lik
                -x(:,nosp)*(S.w_b(nosp).*f_kdt(nosp));
        end
    end %function f_bko


    function [lik dlik]= f_bk(RateParams)                           % get lik and grad

        xk      = RateParams'*Sim.x;                                % filtered stimulus
        hs      = zeroy;                                            % incorporate spike history terms
        for l=1:Sim.M
            hs  = hs+E.omega*S.h(:,:,l);
        end
        s       = repmat(xk,Sim.N,1) + hs;

        f_kdt   = exp(s)*Sim.dt;                                    % shorthand
        lik     = -sum(S.w_b(sp).*log(1-exp(-f_kdt(sp))))...        % liklihood
            +sum(S.w_b(nosp).*f_kdt(nosp));

        if nargout > 1                                              % if gradobj=on
            ef      = exp(f_kdt);                                   % shorthand
            dlik      = x(:,sp)*(S.w_b(sp).*f_kdt(sp)./( ef(sp)-1))... % gradient of lik
                -x(:,nosp)*(S.w_b(nosp).*f_kdt(nosp));
        end
    end %function f_bko

%% MLE for calcium parameters
if Sim.C_params == true
    fprintf('estimating calcium parammeters\n')
    [ve_x fval] = quadprog(M.Q, M.L,[],[],[],[],[0 0 0],[inf inf inf],[1/E.tau_c E.A E.C_0/E.tau_c]+eps);
    Enew.tau_c  = 1/ve_x(1);
    Enew.A      = ve_x(2);
    Enew.C_0    = ve_x(3)/ve_x(1);
    Enew.sigma_c= sqrt(-(0.5*ve_x'*M.Q*ve_x + M.L'*ve_x))/(Sim.T*Sim.dt);
    Enew.lik_c  =fval;
    lik = [lik Enew.lik_c];
end

% % %% MLE for spike history parameters
% % for m=1:Sim.M
% %     Enew.sigma_h(m)= sum(M.v{m})/Sim.T;
% % end
%
% %% MLE for observation parameters
% O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);
% Oind    = find(~isnan(O));
% Enew.sigma_o = sqrt(sum(sum(S.w_b(:,Oind).*(repmat(O(Oind),Sim.N,1)-S.C(:,Oind)).^2)/Sim.T));


%% MLE for observation parameters
if Sim.F_params == true
    fprintf('estimating observation parammeters\n')
    options         = optimset('Display','off','GradObj','off'); %use gradient
    ab_0            = [E.alpha E.beta];
    [ab Enew.lik_o]   = fminunc(@f_ab,ab_0,options);%find MLE
    Enew.alpha=ab(1);
    Enew.beta=ab(2);
    lik = [lik Enew.lik_o];
end

    function lik = f_ab(ab_o)
        Fmean=Hill_v1(E,S.w_b.*S.C);
        A=-[Fmean(:) ones(Sim.N*Sim.T,1)];
        H=A'*A;
        f=A'*repmat(F,Sim.N,1);
        x = quadprog(H,f,[],[],[],[],[0 0],[inf inf]);

    end %function f_ab_0
%%
Enew.lik=sum(lik);
end
