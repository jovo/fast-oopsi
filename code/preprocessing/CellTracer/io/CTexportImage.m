function CTexportImage(proj,frame,which)
filename = proj.filename;
if ~iscell(filename)
    filename = cellstr(filename);
end
outputpath = proj.outputpath;
pathname = proj.pathname;
currentpath = pwd;
cd(outputpath);
    
matfile = [char(filename{frame}),'.mat'];
load(matfile);
viewindex = length(views);
data = views{viewindex};
if which ==1
    im = data.im;
    cd(pathname);
    tiffile = char(filename{frame});
    imwrite(im,tiffile);
else
    mask=data.mask;
    if ~isempty(mask)
        mask = mask(:,:,3);
        mask = uint8(mask) * 255;
        tiffile = [char(filename{frame}),'.seg.tif'];
        imwrite(mask,tiffile);
    end
end
cd(currentpath);          