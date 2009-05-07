%% this file does everything required to compare a few ways of spatial
%% filters
%
% 1) extracts a region of interest in a format useful for matlab
% 2) builds a bunch of filters
% 3) plots what the figures look like
% 3) compares the spike extraction accuracy using the different filters

clear, clc, %close all
foldname    = 'TS108_6';
dir         = '/Users/joshyv/Research/projects/oopsi/meta-oopsi/data/karel/takashi/';
im_dir      = '/Image/';
phys_dir    = '/Physiology/';

for i=12;%:2
    
    % get whole movie
    if i<10                                         % get tif file name
        tifname=[dir foldname im_dir  foldname '_00' num2str(i) '.tif'];
    else
        tifname=[dir foldname im_dir  foldname '_0' num2str(i) '.tif'];
    end
    MovInf  = imfinfo(tifname);                     % get number of frames
    T       = numel(MovInf)/2;                      % only alternate frames have functional data
    DataMat = zeros(MovInf(1).Width*MovInf(1).Height,T);% initialize mat to store movie
    for j=1:2:T*2                                    
        X = imread(tifname,j);%,'PixelRegion'), {ROWS, COLS});
        DataMat(:,(j+1)/2)=X(:);
    end
    
    % manually select ROI
    MeanFrame=mean(DataMat');
    figure(100), clf,
    imagesc(reshape(z1(MeanFrame),MovInf(1).Height,MovInf(1).Width))
    title('select corners, starting top-left, going clock wise')
    [x y]   = ginput(4);
    ROWS    = [round(mean(y(1:2))) round(mean(y(3:4)))];                              % define ROI
    COLS    = [round(mean(x([1 4]))) round(mean(x(2:3)))];                    

    Nrows=diff(ROWS)+1;
    Ncols=diff(COLS)+1;
    ind=zeros(Nrows*Ncols,1);
    j=0;
    for r=1:Nrows
        for c=1:Ncols
            j=j+1;
            ind(j)=MovInf(1).Height*(COLS(1)-1+c)+ROWS(1)+r-1;
        end
    end
    
    % get filters
    R.F=-double(DataMat(ind,:));
    MeanROI=mean(R.F,2);
    figure(101), clf,
    imagesc(reshape(z1(MeanROI),Ncols,Nrows)')
    
%     D = GetSpatialFilter(R.F,Ncols,Nrows);           

    % get spike times
    if i<10
        matname=[dir foldname phys_dir '/T0004E196P1t000' num2str(i) '.mat'];
    else
        matname=[dir foldname phys_dir '/T0006E108P1t00' num2str(i) '.mat'];
    end
    load(matname)
    spt     = Get_spt(Trace.data);                  % get spike times
    R.spt   = round(spt/(length(Trace.data)/T));    % downsample to frame rate

    % plot stuff for consistency check
    sp      = zeros(T,1); sp(R.spt)=1;              % plot spt and F to for consistency check
    figure(7), hold on, stem(sp,'k'), hold off

    dV      = diff(Trace.data);                     % plot spike times and voltage to confirm correctness
    dV      = dV/max(dV);
    sp      = 0*dV;
    sp(spt) = 1;   
    figure(8),  plot(dV), hold on, stem(sp*.5,'k'), hold off
    if i<10                                         % get tif file name
        savname=[foldname '_00' num2str(i) '.mat'];
    else
        savname=[foldname '_0' num2str(i) '.mat'];
    end

    save(['/Users/joshyv/Research/projects/oopsi/spatial-filter/' savname],'D','R')
end

