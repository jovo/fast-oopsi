function Pl = PlotParams(Pl)
%% make plot params

Pl.gray = [.75 .75 .75];        % define gray color
Pl.col  = [1 0 0; 0.2 0.2 1];   % define colors for mean
Pl.ccol = Pl.col+.4;  Pl.ccol(Pl.ccol>1)=1; % define colors for std
Pl.inter= 'none';               % interpreter for axis labels
Pl.fs   = 12;                   % font size
Pl.ms   = 20;                   % marker size for real spike
Pl.vs   = 4;                    % marker size for real spike
Pl.sw   = 1.5;                  % spike width
Pl.lw   = 2;                    % line width

% make xticks
TickSpace = round(Pl.xlims(2)/Pl.nticks);
Pl.XTicks = TickSpace*(1:Pl.nticks);

% labels
Pl.FastLabel = [{'Fast'}; {'Filter'}];
Pl.WienerLabel = [{'Wiener'}; {'Filter'}];
end