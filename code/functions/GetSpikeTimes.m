function spt = GetSpikeTimes(V,thr)
if nargin < 2
    thr = 0.1;
end
window = 5;
% num_samples = 10;

V=V-mean(V);
V=V/max(abs(V));
r=find(V>thr);                          %find time bins during a spike
r1=find(diff(r)>6);                     %find time points marking end of spike
if ~isempty(r1)
    nsp=length(r1);                     %number of spikes - 1
    spt=zeros(nsp,1);                   %initialize spike time vector
    [a b]=max(V(r(1:r1(1))));           %find max of first spike
    spt(1) =  r(b);                     %let first spike time be the time of that max
    for t=2:nsp                         %for each spike
        temp=r(r1(t-1)+1:r1(t));        %time bins during spike
        [a b]=max(V(temp));             %find max
        spt(t) =  temp(b);              %let time be the max
    end
    [a b]=max(V(r(r1(end))+5:r(end)));  %add last spike
    spt = [spt; r(r1(end))+b-1];
else
    [a b]=max(V(r));                   %find max of first spike
    spt(1) =  r(b);                     %let first spike time be the time of that max
end