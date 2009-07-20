function [roc auc] = roc3(x,y)

x   = x/max(abs(x));
ts  = 1:-.1:0;
fn  = zeros(size(ts));
tp  = zeros(size(ts));
fp  = zeros(size(ts));
sp1 = find(y);
sp0 = find(~y);
nsp1 = length(sp1);
nsp0 = length(sp0);
k=0;
for t=ts
    k=k+1;
    tp(k) = sum(x(sp1)>t)/nsp1;
    fp(k) = sum(x(sp1)<t)/nsp1;
    tn(k) = sum(x(sp0)<t)/nsp0;
    fn(k) = sum(x(sp0)>t)/nsp0;
end

roc=tp;
auc=sum(roc);