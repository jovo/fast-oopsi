function roc = roc2(x,y)

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
    roc.tp(k) = sum(x(sp1)>t)/nsp1;
    roc.fp(k) = sum(x(sp1)<t)/nsp1;
    roc.tn(k) = sum(x(sp0)<t)/nsp0;
    roc.fn(k) = sum(x(sp0)>t)/nsp0;
end
