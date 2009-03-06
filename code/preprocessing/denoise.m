clear, clear global, close all
filehead='power20_T';
filefoot='_B.tif';

for n=2:1205
    if n<10
        pad='000';
    elseif n>9 & n<100
        pad='00';
    elseif n>99 & n<1000
        pad='0';
    else
        pad='';
    end
    filename=[filehead pad num2str(n) filefoot];
    folder='C:\power20\';
    total=[folder filename];
    a=imread(total);
    a_crop=a(250:2:449,250:2:449);
    a_cat=a_crop(:);
    Cells(:,n-1)=a_cat;
end

Cells=double(Cells); 
save Cells Cells
MeanFrame=mean(Cells');
MeanMatrix=repmat(MeanFrame,size(Cells,2),1)';
Cells=Cells-MeanMatrix; clear MeanMatrix

siz=size(Cells);

if siz(1)>siz(2)
    SIG=Cells'*Cells;
else
    SIG=Cells*Cells';
end
size(SIG)
[EigTime EigVal]=eig(SIG);

v2=diag(EigVal);
for i=0:length(v2)-1
    SortedEigVals(i+1)=v2(end-i);
end
figure, semilogy(SortedEigVals), title('Eigenvalues')
figure, semilogy(SortedEigVals(1:20)), title('First 20 Eigenvalues')

figure, 
for i=1:20
    subplot(5,5,i), plot(EigTime(:,end-i+1)), title(['Eigentime ' num2str(i)])
    set(gca,'XTick',[]),  axis('tight') %set(gca,'YTick',[]),
end

subplot(6,5,21), plot(EigTime(:,end-49)), title('Eigentime 50'), axis('tight'), set(gca,'XTickLabel',[2 4 6 8 10 12]), 
subplot(6,5,22), plot(EigTime(:,end-99)), title('Eigentime 100'), set(gca,'XTick',[]), axis('tight'), % set(gca,'YTick',[])
subplot(6,5,23), plot(EigTime(:,end-199)), title('Eigentime 200'), set(gca,'XTick',[]), axis('tight'), % set(gca,'YTick',[])
subplot(6,5,24), plot(EigTime(:,end-499)), title('Eigentime 500'), set(gca,'XTick',[]), axis('tight'), % set(gca,'YTick',[])
subplot(6,5,25), plot(EigTime(:,end-999)), title('Eigentime 1000'), set(gca,'XTick',[]), axis('tight'), % set(gca,'YTick',[])

save goods Cells EigTime EigVal v2 MeanFrame
clear all
load goods

root_val=sqrt(EigVal);
xinv=(root_val*EigTime')^-1;
EigSpace=Cells*xinv;

save goods2 EigSpace EigVal EigTime MeanFrame
clear all
load goods2

figure, 
for i=1:20
    subplot(5,5,i), imagesc(reshape(EigSpace(:,end-i+1),100,100))
    colormap('gray'), title(['Eigenspace ' num2str(i)])
    set(gca,'XTick',[]), set(gca,'YTick',[])
end

subplot(6,5,21), imagesc(reshape(EigSpace(:,end-49),100,100)), colormap('gray'), title('Eigentime 50'), set(gca,'XTick',[]), set(gca,'YTick',[])
subplot(6,5,22), imagesc(reshape(EigSpace(:,end-99),100,100)), colormap('gray'), title('Eigentime 100'), set(gca,'XTick',[]), set(gca,'YTick',[])
subplot(6,5,23), imagesc(reshape(EigSpace(:,end-199),100,100)), colormap('gray'), title('Eigentime 200'), set(gca,'XTick',[]), set(gca,'YTick',[])
subplot(6,5,24), imagesc(reshape(EigSpace(:,end-499),100,100)), colormap('gray'), title('Eigentime 500'), set(gca,'XTick',[]), set(gca,'YTick',[])
subplot(6,5,25), imagesc(reshape(EigSpace(:,end-999),100,100)), colormap('gray'), title('Eigentime 1000'), set(gca,'XTick',[]), set(gca,'YTick',[])

DeNoised=EigSpace(:,end-13:end)*EigVal(end-13:end,end-13:end)*EigTime(:,end-13:end)';
save DeNoised DeNoised MeanFrame 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make video
%% no noise 
%% no mean
clear all
load DeNoised

mx=max(max(DeNoised));
mn=min(min(DeNoised));
DeNoised=uint8(floor(1+255*(DeNoised-mn)/(mx-mn)));

cmap=[0:1/256:1-1/256'; 0:1/256:1-1/256'; 0:1/256:1-1/256';]';
for i=1:size(DeNoised,2)
    mov(i)=im2frame(reshape(DeNoised(:,i),100,100),colormap(cmap));
end

filename='NoNoiseNoMean';
aviobj=avifile(filename,'compression','None','fps',15);
aviobj = addframe(aviobj,mov);
aviobj=close(aviobj);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make video
%% yes noise
%% no mean
clear all
load goods

mx=max(max(Cells));
mn=min(min(Cells));
Cells=uint8(floor(1+255*(Cells-mn)/(mx-mn)));

for i=1:size(Cells,2)
    mov(i)=im2frame(reshape(Cells(:,i),100,100),colormap('gray'));
end

filename='NoMean';
aviobj=avifile(filename,'compression','None','fps',15);
aviobj = addframe(aviobj,mov);
aviobj=close(aviobj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make raw video
%% yes noise
%% yes mean
clear all
load Cells

mx=max(max(Cells));
mn=min(min(Cells));
Cells=uint8(floor(1+255*(Cells-mn)/(mx-mn)));

for i=1:size(Cells,2)
    mov(i)=im2frame(reshape(Cells(:,i),100,100),colormap('gray'));
end

filename='Raw';
aviobj=avifile(filename,'compression','None','fps',15);
aviobj = addframe(aviobj,mov);
aviobj=close(aviobj);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make video
%% no noised
%% yes mean
clear all
load DeNoised

mx=max(max(DeNoised));
mn=min(min(DeNoised));
DeNoised=(DeNoised-mn)/(mx-mn);
MeanFrame=(MeanFrame-min(MeanFrame))/(max(MeanFrame)-min(MeanFrame));
DeNoised=DeNoised+repmat(MeanFrame(:),1,size(DeNoised,2));

mx=max(max(DeNoised));
mn=min(min(DeNoised));
DeNoised=uint8(floor(1+255*(DeNoised-mn)/(mx-mn)));

for i=1:size(DeNoised,2)
    mov(i)=im2frame(reshape(DeNoised(:,i),100,100),colormap('gray'));
end

filename='NoNoise';
aviobj=avifile(filename,'compression','None','fps',15);
aviobj = addframe(aviobj,mov);
aviobj=close(aviobj);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make video double-precision
%% no noise 
%% no mean
% clear all
% load DeNoised
% 
% mx=max(max(DeNoised));
% mn=min(min(DeNoised));
% DeNoised=1+(2^16-1)*(DeNoised-mn)/(mx-mn);
% 
% cmap=[0:1/2^16:1'; 0:1/2^16:1'; 0:1/2^16:1';]';
% % cmap=[1:2^16-1'; 1:2^16-1'; 1:2^16-1';]';
% for i=1:size(DeNoised,2)
%     mov(i)=im2frame(reshape(DeNoised(:,i),100,100),colormap(cmap));
% end
% 
% filename='NoNoiseNoMean';
% aviobj=avifile(filename,'compression','None','fps',15);
% aviobj = addframe(aviobj,mov);
% aviobj=close(aviobj);
