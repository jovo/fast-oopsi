function n = SubSampleSpikes(D,freq)
% this function subsamples spikes by the factor 'freq', 
% ie, if freq=2, and there were 500 image frames, 
% then the length of n is 1000
% we assume that no frames have been skipped/dropped.

Tstart  = find(diff(D.TTL)>1000,5);             %find TTL onset times, the 5th onset time is the actual first frame
n       = zeros(D.T_o*freq,1);                     %extize spike train
spt     = (D.spt-Tstart(end))/D.rate;           %find spike times, and subtract ephys data before imaging started

MedFrameDur = mean(diff(D.FrameStartTime));   %find the median of the frame duration
dt = MedFrameDur/freq;                          %make desired dt
FrameStartTime = (0:dt:(D.T_o*freq)*dt)';       %make a fake frame start time

for k=1:length(spt)                             %step thru each interval
        sub=find(spt(k)<FrameStartTime,1);      %add 1 to n at that time step if a spike happened within that interval
        if spt(k)<FrameStartTime(sub)-dt/2, 
            s=sub-1; 
        else s=sub;
        end
        n(s)=n(s)+1;
end

T=min(find((length(D.V)-Tstart(end))/D.rate<FrameStartTime,1),length(n)); %find the frame of the last ephys data
n=n(1:T);                                       %truncate the spike train to only go as far as the ephys wentend