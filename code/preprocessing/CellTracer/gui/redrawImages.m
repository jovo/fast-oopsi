function [left right] = redrawImages(im,mask,ccolor)
currentview = getCurrentView(gcf);
if currentview ~= 3
    msgbox('redrawImages::program bug!');
else
    handles = guidata(gcf);
    frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
    bwim = getCurrentLabels(gcf,frame);
    
    bwmask =  getCurrentLabels(gcf,frame+1);
    if isempty(bwmask) 
        bwmask = mask;
        bwmask = bwlabel(bwmask,4); 
    end
    %redraw the left axis based on view
    ax1 = getAxes(gcf,1); subplot(ax1); v = axis;
    left = drawimage(currentview,1,im,ccolor,bwim);
    imshow(left);
    %redraw the right axis based on view
    ax2 = getAxes(gcf,2); subplot(ax2);
    right = drawimage(currentview,2,mask,ccolor,bwmask);
    imshow(right);
    axis(v);
end

