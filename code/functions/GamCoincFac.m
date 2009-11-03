%
%    G = GamCoincFac(PredSpkTrain,TargetSpkTrain, SamplingFreq)
%       Calculates the Gamma factor for the Target Spike Train
%       TargetSpkTrain and the Modeled Spike Train PredSpkTrain.
%       SamplingFreq the sampling frequency in Hz, and the spike trains are
%       a list of spike indices (spike times = SpkTrain / SamplingFreq).  
%   
%   See
%       Kistler et al, Neural Comp 9:1015-1045 (1997)
%       Jolivet et al, J Neurophysiol 92:959-976 (2004)
%   for further details
%
%
%           - Renaud Jolivet 2007.05.
%           - modified by R. Naud 2007.09. (siingularity and break in loop
%          and empty PredSpkTrain handling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function G = GamCoincFac(PredSpkTrain,TargetSpkTrain,SamplingFreq)

if isempty(PredSpkTrain), G = 0; return, end


% Some parameters
DeltaWindow     =   2e-03;                      % [sec]         % half-width of coincidence detection (2 msec)
DeltaBins       =   DeltaWindow*SamplingFreq;   % [time bins]   % half-width of coincidence detection (2 msec)
NSpikesPred     =   length(PredSpkTrain);
NSpikesTarget   =   length(TargetSpkTrain);
%
% Compute frequencies, normalisation and average coincidences 
FreqPred        =   SamplingFreq*(NSpikesPred-1)/max((PredSpkTrain(NSpikesPred)-PredSpkTrain(1)),1);
NCoincAvg       =   2*DeltaWindow*NSpikesTarget*FreqPred;
NNorm           =   abs(1-2*FreqPred*DeltaWindow);
%
% Compute the gamma coincidence factor
NCoinc          =   0; 
i               =   1;
while i <= NSpikesTarget
    j=1;
    while j <= NSpikesPred
        if abs(PredSpkTrain(j)-TargetSpkTrain(i)) <= DeltaBins
            NCoinc  =   NCoinc+1;
            i       =   i+1;
            if i> NSpikesTarget, break, end
        end
        j=j+1;
    end
    i=i+1;
end

G               =	(NCoinc-NCoincAvg)/(1/2*(NSpikesPred+NSpikesTarget))*1/NNorm;