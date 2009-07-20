function ROCdata=roc(varargin)
% ROC - Receiver Operating Characteristics.
% The ROC graphs are a useful tecnique for organizing classifiers and
% visualizing their performance. ROC graphs are commonly used in medical
% decision making.
% If you have downloaded partest 
% http://www.mathworks.com/matlabcentral/fileexchange/12705
% the routine will compute several data on test performance.
%
% Syntax: roc(x,alpha)
%
% Input: x - This is the data matrix. The first column is the column of the data value;
%            The second column is the column of the tag: unhealthy (1) and
%            healthy (0).
%          alpha - significance level (default 0.05)
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
elseif nu>3
    error('Warning: Max three input data are required')
end
default.values = {[],0.05,1};
default.values(1:nu) = args;
[x alpha verbose] = deal(default.values{:});
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

lx=length(x(x(:,2)==1)); ly=length(x(x(:,2)==0)); %number of subjects
z=sortrows(x,1);
%find unique values in z
labels=unique(z(:,1));
ll=length(labels); %count unique value
a=zeros(ll,2); %array preallocation
xbar=mean(x(x(:,2)==1),1); ybar=mean(x(x(:,2)==0),1);
for K=1:ll
     I=length(z(z(:,1)<=labels(K))); %Set unique value(i) as cut-off
     table(1)=length(z(z(1:I,2)==1)); %true positives
     table(2)=length(z(z(1:I,2)==0)); %false positives
     test=[table;[lx ly]-table]; %complete the table
     if xbar>ybar
         test=fliplr(test); %mirror the ROC curve
     end
     a(K,:)=(diag(test)./sum(test)')'; %Sensitivity and Specificity
end
xroc=1-a(:,2); yroc=a(:,1); %ROC points

Area=1-trapz(yroc,xroc); %estimate the area under the curve
%standard error of area
Area2=Area^2; Q1=Area/(2-Area); Q2=2*Area2/(1+Area);
V=(Area*(1-Area)+(lx-1)*(Q1-Area2)+(ly-1)*(Q2-Area2))/(lx*ly);
Serror=realsqrt(V);
if verbose
    %confidence interval 
    cv=realsqrt(2)*erfcinv(alpha);
    ci=Area+[-1 1].*(cv*Serror);
    %z-test
    SAUC=(Area-0.5)/Serror; %standardized area
    p=1-0.5*erfc(-SAUC/realsqrt(2)); %p-value
    %Performance of the classifier
    if Area==1
        str='Perfect test';
    elseif Area>=0.90 && Area<=0.99
        str='Excellent test';
    elseif Area>=0.80 && Area<=0.89
        str='Good test';
    elseif Area>=0.70 && Area<=0.79
        str='Fair test';
    elseif Area>=0.60 && Area<=0.69
        str='Poor test';
    elseif Area>=0.50 && Area<=0.59
        str='Fail test';
    else
        str='Failed test - less than chance';
    end
    %display results
    disp('ROC CURVE ANALYSIS')
    disp(' ')
    tr=repmat('-',1,80);
    disp(tr)
    str2=['AUC\t\t\t\tS.E.\t\t\t\t\t' num2str((1-alpha)*100) '%% C.I.\t\t\t\tComment\n'];
    fprintf(str2)
    disp(tr)
    fprintf('%0.5f\t\t\t%0.5f\t\t\t%0.5f\t\t\t%0.5f\t\t\t%s\n',Area,Serror,ci,str)
    disp(tr)
    fprintf('Standardized AUC\t\t1-tail p-value\n')
    fprintf('%0.4f\t\t\t\t\t%0.6f',SAUC,p)
    if p<=alpha
        fprintf('\t\tThe area is statistically greater than 0.5\n')
    else
        fprintf('\t\tThe area is not statistically greater than 0.5\n')
    end
    %display graph
    plot([0;xroc],[0;yroc],'r.-')
    hold on
    plot([0 1],[0 1],'k')
    hold off
    xlabel('False positive rate (1-Specificity)')
    ylabel('True positive rate (Sensitivity)')
    title('ROC curve')
    axis square
    %if partest.m was downloaded
    if p<=alpha
        %the best cut-off point is the cloSerrorst point to (0,1)
        d=realsqrt(xroc.^2+(1-yroc).^2); %apply the Pitagora's theorem
        [Y,J]=min(d); %find the least distance
        co=labels(J); %Set the cut-off point
        hold on
        plot(xroc(J),yroc(J),'bo')
        hold off
        disp(' ')
        fprintf('cut-off point for best Sensitivity and Specificity (blu circle in plot)= %0.4f\n',co)
        disp('Press a key to continue'); pause
        %table at cut-off point
        I=length(z(z(:,1)<=co));
        table(1)=length(z(z(1:I,2)==1));
        table(2)=length(z(z(1:I,2)==0));
        test=[table;[lx ly]-table];
        if xbar>ybar
            test=fliplr(test);
        end
        disp('Table at cut-off point')
        disp(test)
        disp(' ')
        try        
            partest(test)
        catch ME
            disp(ME)
            disp('If you want to calculate the test performance at cutoff point please download partest.m from Fex')
            disp('http://www.mathworks.com/matlabcentral/fileexchange/12705')
        end
    end
end
if nargout
    ROCdata.AUC=Area;
    ROCdata.SE=Serror;
    ROCdata.xr=xroc;
    ROCdata.yr=yroc;
end
