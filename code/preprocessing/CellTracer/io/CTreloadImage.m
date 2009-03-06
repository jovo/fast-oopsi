function CTreloadImage(project,h)
    ghandles = guidata(h); 
    filename = project.filename;
    if ~iscell(filename)
        filename = cellstr(filename);
    end
    ax1 = subplot('position',[0.005, 0.20, 0.375 0.9]); 
    ax2 = subplot('position',[0.385, 0.20, 0.375 0.9]);
    try
        img = imread([project.pathname,char(filename{1})]);
        ax1 = getAxes(gcf,1);subplot(ax1);
        imshow(img,[]); 
        ax2 = getAxes(gcf,2);subplot(ax2);
        imshow(uint8(img <0),[]);
        linkaxes([ax1 ax2]);
    catch 
    end
    guidata(h,ghandles);
