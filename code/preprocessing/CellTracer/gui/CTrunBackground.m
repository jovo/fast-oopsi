function CTrunBackground(which)
h = gcf;
handles = guidata(h);
para = CTgetParameters;

handles.undo.flag = 1;
handles.undo.data = handles.im{handles.io.index};
handles.undo.type = 'segmentation';
updateMenuStatus();

if para.allframes > 0
    frames = length(handles.im{handles.io.index});
    for frame = 1: frames
        set(handles.StatusBar,'String', ['Working on frame: [' num2str(frame) ']']);
        set(handles.FrameIndex,'String',num2str(frame));
        showFrame(h,frame);
        [im mask] = getCurrentImage(h,[]);
        temp = im; 
        if para.sequential >0 
            temp(mask(:,:,2) > 0) = 255;
            tempmask = mask(:,:,2);
        else
            tempmask = [];
        end
        switch which
            case {1}
                datamask = CTbackground1(temp,{para.backgroundspread,para.disk}, para.debug);
            case {2}
                datamask = CTbackground2(temp,{para.mincellvolume,para.maxcellvolume, ...
                        para.disk, para.which}, para.debug);
            case {3}
                datamask = CTbackground3(temp,{para.maxhalfwidth,para.backgroundspread, ...
                        para.disk, para.greedy,para.fillholes}, para.debug);
            case {4}
                datamask = CTbackground4(temp,mask(:,:,2),{para.maxhalfwidth,para.backgroundspread, ...
                        para.minhalfwidth}, para.debug);
            otherwise
        end
        
        
        if para.sequential >0 
            datamask(mask(:,:,2) > 0) = 0;
        end
        mask(:,:,2) = datamask==0;
        setCurrentImage(h,im,mask,frame);
        figure(handles.figure1);
        showFrame(h,[]);
    end
    set(handles.StatusBar,'String', 'All Done');
else
    set(handles.StatusBar,'String', 'Working on current frame');
    [im mask] = getCurrentImage(h,[]);
    temp = im; 
    if para.sequential >0 
        temp(mask(:,:,2) > 0) = 255;
        tempmask = mask(:,:,2);
    else
        tempmask = [];
    end
    switch which
        case {1}
            datamask = CTbackground1(temp,{para.backgroundspread,para.disk}, para.debug);
        case {2}
            datamask = CTbackground2(temp,{para.mincellvolume,para.maxcellvolume, ...
                    para.disk, para.which}, para.debug);
        case {3}
            datamask = CTbackground3(temp,{para.maxhalfwidth,para.backgroundspread, ...
                        para.disk, para.greedy,para.fillholes}, para.debug);
        case {4}
            datamask = CTbackground4(temp,tempmask,{para.maxhalfwidth,para.backgroundspread, ...
                        para.minhalfwidth}, para.debug);
        otherwise
    end
    
    if para.sequential >0 
        datamask(mask(:,:,2) > 0) = 0;
    end
    if para.debug >0 
        rgb(:,:,1) = im;
        rgb(:,:,2) = uint8(datamask==0)*255;
        rgb(:,:,3) = 0;
        figure,imshow(rgb);
    else
        mask(:,:,2) = datamask==0;
        setCurrentImage(h,im,mask,[]);
        figure(handles.figure1);
        showFrame(h,[]);
    end
    set(handles.StatusBar,'String', 'Done');  
end