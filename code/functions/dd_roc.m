function [e,thr] = dd_roc(w,a,b,frac_rej)
% e = dd_roc(W,A,B,frac_rej)
%
% Find for a (data description) method W (trained with A) the
% Receiver Operating Characteristic curve over dataset B. The user can
% supply a vector of target rejection rates for which the ROC curve
% can be computed. These target rejection rates are then obtained for
% the training set A, and applied to the testing set B.
%
% The results are returned in e.  The first column gives the fraction
% of target objects rejected, the second column the fraction of
% outlier objects accepted.
%
% e = dd_roc(W,A)
%
% When no seperate training and testing set are given, the threshold
% is just varied over A, and the target rejection and outlier
% acceptance rates are calculated.
%
% see also dd_auc, dd_error.

% Copyright: D. Tax, R.P.W. Duin, davidt@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

if (nargin<4)
  frac_rej = [0.01 0.025 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.5];
end
if (nargin==3) & (~isa(b,'dataset'))
  frac_rej = b;
  b = a;
end
if (nargin<3)
  b = a;
  frac_rej = 50;
end


% if frac_rej>1 then we just want to change the threshold, and we do
% not require a specific target rejection rate:
if max(length(frac_rej))==1 & frac_rej>1
  nrfrac = frac_rej;
  frac_rej = [];
else
  nrfrac = length(frac_rej);
end
e = zeros(nrfrac,2);

% all checks
if isa(w,'dataset')
  error('Please make W a mapping');
end

% map the labels of a and b to faster numerical labels:
%[nlaba,nlabb,lablist] = renumlab(getlab(a),getlab(b));
%a = dataset(+a,nlaba);
%b = dataset(+b,nlabb);

[W,lablist,map,k,c,v,par] = mapping(w);

if (isempty(k) | (k*c) == 0)  % the mapping is untrained

  if strcmp(map,'svdd') | strcmp(map,'newsvdd')
    % This is an exception, because SVDD has to
                           % be trained for each threshold.
    sigm = range_svdd(a,frac_rej,W{2});
    [I1,I2] = find_target(b);                       
      
      w=newsvdd(a,[],sigm(1),W{2});
    %    disp('Train the SVDD over and over again...');
    for i=1:nrfrac
      w = newsvdd(a,[],sigm(i),W{2});
      d = +(w*b);
      d = d(:,1)-d(:,2); % correct for threshold
      e(i,1) = sum(d(I1,1)<0);   % target objects rejected
      e(i,2) = sum(d(I2,1)>=0);    % outlier objects accepted
    end
    e(:,1) = e(:,1)/length(I1);
    e(:,2) = e(:,2)/length(I2);
    return                 % we are done here.
  else
    if strcmp(map,'newsvdd1')
  disp('newsvdd1: adjusting \nu')
      [I1,I2] = find_target(b);                       
        
      %    disp('Train the SVDD over and over again...');
      for i=1:nrfrac
        w = newsvdd1(a,W{1},frac_rej(i));
        d = +(w*b);
        d = d(:,1)-d(:,2); % correct for threshold
        e(i,1) = sum(d(I1,1)<0);   % target objects rejected
        e(i,2) = sum(d(I2,1)>=0);    % outlier objects accepted
      end
      e(:,1) = e(:,1)/length(I1);
      e(:,2) = e(:,2)/length(I2);
      return                 % we are done here.
    else
      w = w*a;               % the mapping has to be trained
    end
  end
end

% w should be a trained mapping, now compute the thresholds:
if strcmp(map,'nndd')    % Another exception, use leave-one-out from NNDD
  [W,lablist,map,k,c,v,par] = mapping(w);
  d =  dist2dens(log(W.fit),log(W.scale));
else
  d = w*a;               %map the training set
end

% the first class if per definition the target class:
d = sort(+d(:,1));
if (d(end)<1e-10)
  dolog = 1;
  d = log(d);
else
  dolog=0;
end
% now, do we want specific target rejection rates, or just change the
% threshold?
if isempty(frac_rej)
  thr = linspace(d(1),d(end),nrfrac);
else
  frac = round(frac_rej*length(d));
  thr = zeros(length(nrfrac),1);  % computation thresholds:
  for i=1:nrfrac
    if (frac(i)==0)
      thr(i) = d(1);
    else
      thr(i) = (d(frac(i)) + d(frac(i)+1))/2;
    end
  end
end

% use the threshold to test on set b:
[I1,I2] = find_target(b);
d = w*b;
if dolog
  d = log(d);
end
for i=1:nrfrac
  e(i,1) = sum(d(I1,1)<thr(i));    % target objects rejected
  e(i,2) = sum(d(I2,1)>=thr(i));   % outlier objects accepted
end
e(:,1) = e(:,1)/length(I1);
e(:,2) = e(:,2)/length(I2);

return
