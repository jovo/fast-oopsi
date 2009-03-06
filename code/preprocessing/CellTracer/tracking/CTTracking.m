function [forwardscores,backwardscores] = CTTracking(labels1, labels2,alignment, para)
n1 = para{1};n2 = para{2};n3 = para{3};t1 = para{4};

[r1 c1] = find(labels1>0);
rmax= max(r1); rmin = min(r1);
cmax= max(c1); cmin = min(c1);
a = labels1; a = a(rmin:rmax,cmin:cmax);
neighbours{1} = GetNeighbours(a,n1,n2,3);

[r2 c2] = find(labels2>0);
rmax= max(r2); rmin = min(r2);
cmax= max(c2); cmin = min(c2);
a = labels2; a = a(rmin:rmax,cmin:cmax);
neighbours{2} = GetNeighbours(a,n1,n2,3);

rmax= max([max(r1) max(r2)]); rmin = min([min(r1) min(r2)]);
cmax= max([max(c1) max(c2)]); cmin = min([min(c1) min(c2)]);
a = labels1; a = a(rmin:rmax,cmin:cmax);
b = labels2; b = b(rmin:rmax,cmin:cmax);
if ~isempty(alignment)
    [a b] = CTshowGlobalAlignment(a,b,alignment,0);
end
method = -3; %back compatibility, all dataanalysis prior to March 21st, 2008 use default value -3

links = CTcorrespondence(method,a,b,neighbours{1},neighbours{2},n3,1,t1,0);
[newf newb] = CTresolveCorrespondence(links(:,:,1),links(:,:,2)');
forwardscores = newf;
backwardscores = newb;

