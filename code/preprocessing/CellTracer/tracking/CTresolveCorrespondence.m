function [newf newb] = CTresolveCorrespondence(forward, backward)
f1 = forward;
b1 = backward;
common = f1 > 0 & b1' > 0;
f1(common == 0) = 0;
b1(common' == 0) = 0;

%find the points that are both row maximum and column maximums;
%repeat this after removing these points until no more points can be
%removed
notdone = 1;
[m n] = size(f1);
one_mask = zeros(m,n,'uint8');
while notdone >0
    data = f1 .* b1'; %(f1+b1').*(f1 +b1');
    rmax = repmat(max(data,[],1),m,1);
    cmax = repmat(max(data,[],2),1,n);
    f1(data >= rmax & data >= cmax & rmax >0) = 1;
    b1(data' >= rmax' & data' >= cmax' & rmax' > 0) = 1;
    data = f1 .* b1';

    rmax = max(data,[],1);
    mask = data >=1;
    one_index = find(rmax == 1);
    f1(:,one_index) = 0; 
    b1(one_index,:) = 0; 
    one_mask(mask) = 1;
    notdone = sum(mask(:));
end
newf = forward; newf(one_mask < 1) = 0;
newb = backward; newb(one_mask' < 1) = 0;
