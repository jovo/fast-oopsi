function err = dd_auc(e,bnd)
% DD_AUC compute the integrated error under the ROC curve
%
%   err = dd_auc(e,bnd)
%
% Compute the integrated error under the ROC curve, given in e.
% The user can supply a lower and upper bound over which the ROC is integrated:
% bnd = [lowerbnd, upperbnd], for instance bnd = [0.05 0.5].
%
% see also dd_error, dd_roc.

% Copyright: D. Tax, R.P.W. Duin, davidt@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

if (nargin<2)
  bnd = [0.05 0.5];
end

% make e the correct size:
if (size(e,2)~=2)
  e = e';
  if (size(e,2)~=2)
    error('Please make e an Nx2 matrix');
  end
end

% extend the e to the bounds if necessary:
% first, is there something smaller than the lower bound?
I = find(e(:,1)<bnd(1));
% extend to low target rejection rates:
if isempty(I)
  e = [bnd(1) 1; min(e(:,1)) 1; e];
else
  % else interpolate:
  if (length(I)==size(e,1))
    e = [bnd(1) 1];
  else
    x1 = e(I(end),1); y1 = e(I(end),2);
    x2 = e(I(end)+1,1); y2 = e(I(end)+1,2);
    y = y1 + (y2-y1)*(bnd(1)-x1)/(x2-x1);
    e(I,:) = [];
    e = [bnd(1) y; e];
  end
end
I = find(e(:,1)>bnd(2));
if isempty(I)   % extend to high target rejection rates
  e = [e; bnd(2) min(e(:,2))];
else                     % interpolate:
  if (length(I)==size(e,1))
    e = [bnd(1) min(e(:,2))];
  else
    x1 = e(I(1)-1,1); y1 = e(I(1)-1,2);
    x2 = e(I(1),1); y2 = e(I(1),2);
    y = y1 + (y2-y1)*(bnd(2)-x1)/(x2-x1);
    e(I,:) = [];
    e = [e; bnd(2) y];
  end
end

% do we have to enforce that the e(:,1) and e(:,2) are strictly
% increasing/decreasing??

%hold on;plot(e(:,1),e(:,2),'g-');
% integrate:
err = 0;
de = diff(e);
for i=1:size(e,1)-1
  err = err + de(i,1)*(e(i+1,2) - de(i,2)/2);
end

return
