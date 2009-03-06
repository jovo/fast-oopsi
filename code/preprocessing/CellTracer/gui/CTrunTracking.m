function CTrunTracking(which)
h = gcf;
handles = guidata(h);
para = CTgetParameters;

handles.undo.flag = 1;
updateMenuStatus();

if para.allframes > 0
    frames = handles.io.frames;
    for frame = 1: frames-1
        set(handles.StatusBar,'String', ['Working on frame: [' num2str(frame) ']']);
        pause(0.1);
        set(handles.FrameIndex,'String',num2str(frame));
        showFrame(h,frame);
        [im mask] = getCurrentImage(h,[]);
        switch which
            case {1}
                alignment = CTglobalAlignment(im>0, mask>0,{para.shift},0);
                if para.debug > 0
                    CTshowGlobalAlignment(im>0,mask>0,alignment,100);
                else
                    [clink ccolor oldalignment] = getCurrentLink(h,frame);
                    setCurrentLink(h,clink,ccolor,alignment,frame);
                end
            case {2}
                frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
                labels1 = getCurrentLabels(h,frame);
                labels2 = getCurrentLabels(h,frame+1);
                if isempty(labels1) || isempty(labels2)
                    set(handles.StatusBar,'String', 'Cannot perform cell tracking operation.');  
                else
                    [clink ccolor alignment] = getCurrentLink(h,frame);  %get all verified links
                    [forwardscores,backwardscores] = CTTracking(labels1,labels2,alignment, {para.shift,...
                        para.shiftscalefactor, para.celldisplacement,para.minscore});
                    clink = backwardscores > 0;
                    clink = clink';

                    ccolor = CTassignLinkColor(clink, []);
                    setCurrentLink(h,clink,ccolor,alignment,frame);
                    showFrame(h,frame);
                end
            otherwise
        end
    end
    set(handles.StatusBar,'String', 'All Done');
else
    set(handles.StatusBar,'String', 'Working on current frame');
    [im mask] = getCurrentImage(h,[]);
    switch which
        case {1}
            alignment = CTglobalAlignment(im>0, mask>0,{para.shift},0);
            if para.debug > 0
                CTshowGlobalAlignment(im>0,mask>0,alignment,100);
            else
                [clink ccolor oldalignment] = getCurrentLink(h,[]);
                setCurrentLink(h,clink,ccolor,alignment,[]);
                CTshowGlobalAlignment(im>0,mask>0,alignment,100);
            end
        case {2}
            frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
            labels1 = getCurrentLabels(h,frame);
            labels2 = getCurrentLabels(h,frame+1);
            if isempty(labels1) || isempty(labels2)
                set(handles.StatusBar,'String', 'Cannot perform cell tracking operation.');  
            else
                [clink ccolor alignment] = getCurrentLink(h,frame);  %get all verified links
                [forwardscores,backwardscores] = CTTracking(labels1,labels2,alignment, {para.shift,...
                    para.shiftscalefactor, para.celldisplacement,para.minscore});
                clink = backwardscores > 0;
                clink = clink';
                
                ccolor = CTassignLinkColor(clink, []);
                setCurrentLink(h,clink,ccolor,alignment,frame);
                showFrame(h,frame);
            end
             
        otherwise
    end
    set(handles.StatusBar,'String', 'Done');  
end