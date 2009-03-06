function [numframes statustext] = getImages(filename,pathname,outputpath,h)
    ghandles = guidata(h); 
    if iscell(filename)
        frames  = length(filename);
    else
        frames = 1;
        filename = cellstr(filename);
    end
    channel = 0;
    ax1 = subplot('position',[0.005, 0.20, 0.375 0.9]); 
    ax2 = subplot('position',[0.385, 0.20, 0.375 0.9]);
    try
        img = imread([pathname,char(filename{1})]);
        if length(size(img)) == 3
            ax1 = getAxes(gcf,1);subplot(ax1);
            imshow(img,[]); 
            channel = selectChannel();
        end
    catch 
        im = {};        
        statustext = ['Reading image files ... failed']; 
    end
    if ~exist(outputpath,'dir')
        mkdir(outputpath);
    end
    mmin = []; mmax = [];
    currentpath = pwd;
    try
        for frame = 1: frames
            img = imread([pathname,char(filename{frame})]);
            if channel == 0 %single channel
                %do nothing
            elseif channel >3 %rgb to grey
                img = rgb2gray(img);
            else
                img = img(:,:,channel);
            end
            if isempty(mmin)
                mmin = min(img(:)); mmax = max(img(:));
            else
                mmin = min(mmin,min(img(:)));
                mmax = max(mmax,max(img(:)));
            end
            ax1 = getAxes(gcf,1);subplot(ax1);
            imshow(img,[]); pause(0.1);
            ax2 = getAxes(gcf, 2);subplot(ax2);
            temp = uint8(img(:,:,1)<0)*255;
            imshow(temp);
            linkaxes([ax1 ax2]);
            data.source = char(filename{frame});
            data.tag = 'input';
            data.im = img;
            data.mask = [];
            cd(outputpath);
            matfile = [data.source,'.mat'];
            views{1} = data;
            save(matfile, 'views');
            cd(currentpath);
        end
        cd(outputpath);
        for frame = 1: frames
            matfile = [char(filename{frame}),'.mat'];
            load(matfile);
            data = views{1};
            temp = uint8((double(data.im) - double(mmin)) / double(mmax - mmin) * 255.0); temp(temp==0) = 1;
            data.im = temp;
            views{1} = data;
            save(matfile, 'views');
            %im{frame}.raw = temp;
        end
        cd(currentpath);
        statustext = ['Reading image files ... done. Total number of images: ', num2str(frames)]; 
        numframes = frames;
    catch 
        %im = {};
        numframes = 0;
        cd(currentpath);
        filename(frame,:);
        statustext = ['Reading image files ... failed to read file: ',filename{frame}]; 
    end
    guidata(h,ghandles);
function result = selectChannel()
    str = {'Red','Green','Blue','Gray'};
    [s,v] = listdlg('PromptString','Select a channel:',...
                'SelectionMode','single',...
                'ListString',str,...
                'InitialValue',1,...
                'Name', 'Channel Selection',...
                'ListSize',[200 80]);
    if isempty(s)
        s = 1;
    end
    result = s;     
