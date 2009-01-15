function spt = Get_spt(V)
% get spike times sampled at ephys sampling rate

if any(V>0)                             %if there is a spike
    [r c v]=find(V>0);                  %find time bins during a spike
    [r1 c1 v1]=find(diff(r)>6);         %find time points marking end of spike
    nsp=length(r1);                     %number of spikes - 1
    spt=zeros(nsp,1);                   %initialize spike time vector
    [a b]=max(V(r(1:r1(1))));           %find max of first spike
    spt(1) =  r(b);                     %let first spike time be the time of that max
    for t=2:nsp                         %for each spike
        temp=r(r1(t-1)+1:r1(t));        %time bins during spike
        [a b]=max(V(temp));             %find max
        spt(t) =  temp(b);              %let time be the max
    end
    [a b]=max(V(r(r1(end)):r(end)));    %add last spike
    spt = [spt; r(r1(end))+b-1];
else
    spt=[];
end