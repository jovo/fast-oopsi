function CTrunInitialPartition(which,type)
h = gcf;
handles = guidata(h);
para = CTgetParameters;

handles.undo.flag = 1;
%handles.undo.data = handles.im{handles.io.index};
%handles.undo.type = 'segmentation';
updateMenuStatus();

if para.allframes > 0
    frames = handles.io.frames;
    for frame = 1: frames
        set(handles.StatusBar,'String', ['Working on frame: [' num2str(frame) ']']);
        set(handles.FrameIndex,'String',num2str(frame));
        showFrame(h,frame);
        [im mask] = getCurrentImage(h,[]);
        if para.sequential == 0 
            mask(mask>0) = 0;
        end
        switch which
            case {1}
                datamask = CTbackground1(im,mask,{para.backgroundspread,para.disk,para.padding}, para.debug);
            case {2}
                datamask = CTbackground2(im,mask,{para.mincellvolume,para.maxcellvolume, ...
                        para.disk, para.minintensity,para.maxintensity,para.maxseed}, para.debug);
            case {3}
                datamask = CTbackground3(im,mask,{para.maxhalfwidth,para.backgroundspread, ...
                        para.disk, para.greedy,para.fillholes}, para.debug);
            case {4}
                datamask = CTborder1(im,mask,{para.borderspread,para.minsize,para.borderthreshold,para.maxhalfwidth}, para.debug);
            case {5}
                datamask = CTborder2(im,mask,{para.maxhalfwidth,para.borderspread, ...
                        para.disk, para.rankingmethod,...
                        para.borderthreshold}, para.debug);
            case {6}
                datamask = CTborder3(im,mask,{para.maxhalfwidth,para.borderspread, ...
                        para.disk, para.minvolume,para.rankingmethod,...
                        para.borderthreshold}, para.debug);
            case {7}
                datamask = CTborder4(im,mask,{para.maxsize,para.minpercentile, ...
                        para.disk, para.maxpercentile,para.globalthreshold}, para.debug);
            case {10}
                datamask = CTiterativeSelectiveSegmentation(im,mask,{para.cellscore,para.minvolume, ...
                        para.disk, para.dispercentile}, para.debug);
            case {11}
                datamask = CTiterativeSelectiveSegmentation2(im,mask,{para.maxsize,para.percentile, ...
                        para.disk, para.cellscore,para.celldeg}, para.debug);
            otherwise
        end
        mask(:,:,type) = datamask;
        setCurrentImage(h,im,mask,frame);
        figure(handles.figure1);
        showFrame(h,frame);
    end
    set(handles.StatusBar,'String', 'All Done');
else
    set(handles.StatusBar,'String', 'Working on current frame');
    [im mask] = getCurrentImage(h,[]);
    if para.sequential == 0 
        mask(mask>0) = 0;
    end
    switch which
        case {1}
            datamask = CTbackground1(im,mask,{para.backgroundspread,para.disk,para.padding}, para.debug);
        case {2}
            datamask = CTbackground2(im,mask,{para.mincellvolume,para.maxcellvolume, ...
                    para.disk, para.minintensity,para.maxintensity,para.maxseed}, para.debug);
        case {3}
            datamask = CTbackground3(im,mask,{para.maxhalfwidth,para.backgroundspread, ...
                        para.disk, para.greedy,para.fillholes}, para.debug);
        case {4}
                datamask = CTborder1(im,mask,{para.borderspread,para.minsize,para.borderthreshold,para.maxhalfwidth}, para.debug);
        case {5}
                datamask = CTborder2(im,mask,{para.maxhalfwidth,para.borderspread, ...
                        para.disk, para.rankingmethod,...
                        para.borderthreshold}, para.debug);
        case {6}
                datamask = CTborder3(im,mask,{para.maxhalfwidth,para.borderspread, ...
                        para.disk, para.minvolume,para.rankingmethod,...
                        para.borderthreshold}, para.debug);
        case {7}
                datamask = CTborder4(im,mask,{para.maxsize,para.minpercentile, ...
                        para.disk, para.maxpercentile,para.globalthreshold}, para.debug);
        case {10}
                datamask = CTiterativeSelectiveSegmentation(im,mask,{para.cellscore,para.minvolume, ...
                        para.disk, para.dispercentile}, para.debug);
        case {11}
                datamask = CTiterativeSelectiveSegmentation2(im,mask,{para.maxsize,para.percentile, ...
                        para.disk, para.cellscore,para.celldeg}, para.debug);
        otherwise
    end
    
    if para.debug >0 
        mask(:,:,type) = datamask;
        figure
        imshow(im); hold on
        himage = imshow(uint8(mask)*255);
        set(himage, 'AlphaData', 0.5);
        hold off
    else
        mask(:,:,type) = datamask;
        setCurrentImage(h,im,mask,[]);
        figure(handles.figure1);
        showFrame(h,[]);
    end
    set(handles.StatusBar,'String', 'Done');  
end