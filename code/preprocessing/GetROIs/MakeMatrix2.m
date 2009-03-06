clear, clear global, close all, clc
folder='C:\D\Classes\CSHL\data\Balazs Project#57.OneFolder\';
head='Balazs Project#57_T';
foot='_B.tif';
nFrames=26123;

try
    load ROI_57_1_26123
catch
    fname = GetFname(folder, head, 2, foot);
    ROIdata=imread(fname);
    C = GetROI(ROIdata);
    Cells=zeros(C.nCells,nFrames);
    save ROI_57_1_26123
end

for t=2:nFrames
    try
    fname = GetFname(folder, head, t, foot);
    data=imread(fname);
    
    for n=1:C.nCells
        Cells(n,t-1)=double(sum(data(C.CellPixels{n})))/C.nPixels{n};
    end

    if mod(t,200)==0
        save Data_57_1_26123
        t
    end
    catch, end
end

imagesc(Cells)