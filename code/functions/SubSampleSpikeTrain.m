function n = SubSampleSpikeTrain(frame_onset_times,spt)

n = 0*frame_onset_times;                     %extize spike train
for k=1:length(spt)
    if spt(k)>frame_onset_times(1)
        sub=find(spt(k)<frame_onset_times,1);      %add 1 to n at that time step if a spike happened within that interval
        n(sub)=n(sub)+1;
    end
end