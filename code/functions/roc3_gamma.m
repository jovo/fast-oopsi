function ROCdata=roc3_gamma(varargin)
% ROC - Receiver Operating Characteristics.
% The ROC graphs are a useful tecnique for organizing classifiers and
% visualizing their performance. ROC graphs are commonly used in medical
% decision making.
% If you have downloaded partest
% http://www.mathworks.com/matlabcentral/fileexchange/12705
% the routine will compute several data on test performance.
%
% Syntax: roc(x,alpha,jitter)
%
% Input: x - This is the data matrix. The first column is the column of the data value;
%            The second column is the column of the tag: unhealthy (1) and
%            healthy (0).
%          alpha - significance level (default 0.05)
%          jitter - allowable jitter in bins
%          window - jitter is defined in allowed spike time jitter (if this
%          is true) or as the number of bins to shift the spike times by
%          (if this is false)
%
% Output: The ROC plot;
%         The Area under the curve with Standard error and Confidence
%         interval and comment.
%         Cut-off point for best sensitivity and specificity.
%         (Optional) the test performances at cut-off point.
%
% Example:
%           load rocdata
%           roc(x)
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) ROC curve: compute a Receiver Operating Characteristics curve.
% http://www.mathworks.com/matlabcentral/fileexchange/19950

%Input Error handling
args=cell(varargin);
nu=numel(args);
if isempty(nu)
    error('Warning: almost the data matrix is required')
elseif nu>5
    error('Warning: Max five input data are required')
end
default.values = {[],0.05,0,1,1};
default.values(1:nu) = args;
[x alpha jitter window verbose] = deal(default.values{:});
if isvector(x)
    error('Warning: X must be a matrix')
end
if ~all(isfinite(x(:))) || ~all(isnumeric(x(:)))
    error('Warning: all X values must be numeric and finite')
end
x(:,2)=logical(x(:,2));
if all(x(:,2)==0)
    error('Warning: there are only healthy subjects!')
end
if all(x(:,2)==1)
    error('Warning: there are only unhealthy subjects!')
end
if ~exist('jitter','var')
    jitter = 0;
end
if nu>=2
    if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
        error('Warning: it is required a numeric, finite and scalar ALPHA value.');
    end
    if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end
    if nu==3
        verbose=logical(verbose);
    end
end
clear args default nu

tr=repmat('-',1,80);
lu=length(x(x(:,2)==1)); %number of unhealthy subjects
lh=length(x(x(:,2)==0)); %number of healthy subjects
z=sortrows(x,1);
%find unique values in z
labels=unique(z(:,1));
ll=length(labels); %count unique value
a=zeros(ll,2); %array preallocation
ubar=mean(x(x(:,2)==1),1); %unhealthy mean value
hbar=mean(x(x(:,2)==0),1); %healthy mean value

% for each spike, make neighboring k bins also spikes
kN = jitter;
if ~exist('window','var') || window == true
    % accept a window of jitters
    sI = find(x(:,2) == 1);
    for K=1:length(sI)
        rng = sI(K)-kN:sI(K)+kN;
        rng = intersect(rng,1:length(x(:,2)));
        x(rng,2) = 1;
    end
else
    % jitter spike times
    x(:,2) = circshift(x(:,2),kN);
end

for K=1:ll
    if hbar<ubar
        % number of times where spike and val > threshold
        TP=length(x(x(:,2)==1 & x(:,1)>labels(K)));
        % number of times no spike and val > threshold
        FP=length(x(x(:,2)==0 & x(:,1)>labels(K)));
        % number of times where spike and val <= threshold
        FN=length(x(x(:,2)==1 & x(:,1)<=labels(K)));
        % number of times whre no spike and val <= threshold
        TN=length(x(x(:,2)==0 & x(:,1)<=labels(K)));
    else
        TP=length(x(x(:,2)==1 & x(:,1)<labels(K)));
        FP=length(x(x(:,2)==0 & x(:,1)<labels(K)));
        FN=length(x(x(:,2)==1 & x(:,1)>=labels(K)));
        TN=length(x(x(:,2)==0 & x(:,1)>=labels(K)));
    end
    a(K,:)=[TP/(TP+FN) TN/(TN+FP)]; %Sensitivity and Specificity
end
xroc=[1; 1-a(:,2); 0]; yroc=[1; a(:,1); 0]; %ROC points

if ~issorted(yroc)
    labels=flipud(labels);
    xroc=flipud(xroc);
    yroc=flipud(yroc);
end

Area=1-trapz(yroc,xroc); %estimate the area under the curve
%standard error of area
Area2=Area^2; Q1=Area/(2-Area); Q2=2*Area2/(1+Area);
V=(Area*(1-Area)+(lu-1)*(Q1-Area2)+(lh-1)*(Q2-Area2))/(lu*lh);
Serror=realsqrt(V);
if verbose
    %confidence interval
    cv=realsqrt(2)*erfcinv(alpha);
    ci=Area+[-1 1].*(cv*Serror);
    if ci(2)>1
        ci(2)=1;
    end
    %z-test
    SAUC=(Area-0.5)/Serror; %standardized area
    p=1-0.5*erfc(-SAUC/realsqrt(2)); %p-value
    %Performance of the classifier
    if Area==1
        str='Perfect test';
    elseif Area>=0.90 && Area<1
        str='Excellent test';
    elseif Area>=0.80 && Area<0.90
        str='Good test';
    elseif Area>=0.70 && Area<0.80
        str='Fair test';
    elseif Area>=0.60 && Area<0.70
        str='Poor test';
    elseif Area>=0.50 && Area<0.60
        str='Fail test';
    else
        str='Failed test - less than chance';
    end
    %display results
    disp(tr)
    fprintf('%0.5f\t\t\t%0.5f\t\t\t%0.5f\t\t%0.5f\t\t\t%s\n',Area,Serror,ci,str)
end
if nargout
    ROCdata.AUC=Area;
    ROCdata.SE=Serror;
    ROCdata.xr=xroc;
    ROCdata.yr=yroc;
end
