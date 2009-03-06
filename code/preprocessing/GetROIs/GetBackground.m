clear, clear global, close all, clc
folder='C:\D\Classes\CSHL\data\Balazs Project#57.OneFolder\';
head='Balazs Project#57_T';
foot='_B.tif';
nFrames=17055;

for t=2:nFrames
    try
    fname = GetFname(folder, head, t, foot);
    data=imread(fname);
    
    BackgroundROI=data(1:100,end-100:end);
    background(t-1)=mean(BackgroundROI(:));

    if mod(t,1000)==0
        save Background_57_1_26123
        t
    end
    catch, end
end