function Out = GetROI(In)
BWs=edge(In, 'prewitt');
se90=strel('line', 3, 90);
se0=strel('line', 3, 0);
BWsdil=imdilate(BWs, [se90 se0]);
BWdfill=imfill(BWsdil,'holes');

[labeled, nObjects] = bwlabel(BWdfill);

nCells=0;
for n=1:nObjects
    pixels=find(labeled==n);
    len=length(pixels);
    if len<20
        labeled(pixels)=0;
    else
        nCells=nCells+1;
        labeled(pixels)=nCells; 
        CellPixels{nCells}=pixels;
        nPixels{nCells}=len;
    end
end

Out.nCells=nCells;
Out.labeled=labeled;
Out.CellPixels=CellPixels;
Out.nPixels=nPixels;