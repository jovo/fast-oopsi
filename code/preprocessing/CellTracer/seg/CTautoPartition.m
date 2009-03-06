function CTautoPartition()
h = gcf;
handles = guidata(h);
%frame= get(handles.FrameIndex,'String'); frame = str2num(frame);
[im mask] = getCurrentImage(h,[]);

if 0> 1
    minsize = 50; maxsize = minsize * 1000; disk = 1; conn = 4;limit = [30 200];
    if isempty(mask)
        %msgbox('Manually select some objects first');
    else
        cells = mask(:,:,3);
        temp = CTiterativeApply(double(cells),0, 4,{@(x) double(x>0).* length(x(x>0))}, []);
        temp = unique(temp(:)); temp = temp(temp>0);
        if isempty(temp)
             %msgbox('Manually select some cells first');

        else
            minsize = temp(1) * 0.5;
            maxsize = temp(end) * 10;
        end

    end
else
    para = CTgetParameters;
    maxsize = para.maxcellvolume;
    minsize = para.mincellvolume;
    disk = para.diskradius;
    conn = 4;
    limit(1) = para.minintensity;
    limit(2) = para.cellentranceintensity;
    limit(3) = para.cellexitintensity;
    limit(4) = para.maxintensity;
end
updateMenuStatus();
a = CTScreen(im,maxsize,minsize,disk,conn,limit,[],[],[]);
guidata(gcf,handles);