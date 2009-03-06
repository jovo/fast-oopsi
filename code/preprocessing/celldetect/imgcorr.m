function [cell_corr,cr,fr]=imgcorr(Img_in,cr,fr,maxchecksize);

%check input arguments

if ~exist('maxchecksize') || prod(size(maxchecksize))==0
maxchecksize=10; %defaults to 10
end

if ~exist('fr') || prod(size(fr))==0
    fr=1; %defaults to 1
end

if ~exist('cr') || prod(size(cr))==0  %if cr is non-existent,cr is calculated below
    
        stdev=zeros(1,maxchecksize);
            for cr=1:maxchecksize
                h_ext=zeros(2*cr+1+2*fr,2*cr+1+2*fr);
                h_ext((fr+1):(2*cr+1+fr),(fr+1):(2*cr+1+fr))=fspecial('disk',cr);
                C = normxcorr2(h_ext,real(Img_in));
                C(C<0)=0;
                stdev(cr)=std2(C);
            end

        [maxval,maxind]=max(stdev);
        cr=maxind;
end
canvas=imresize(real(Img_in),size(Img_in)+2*(cr+fr));
canvas((cr+fr+1):((size(Img_in,1)+cr+fr)),(cr+fr+1):((size(Img_in,2)+cr+fr)))=real(Img_in);
h_ext=zeros(2*(cr+fr)+1,2*(cr+fr)+1);
h_ext(fr+(1:(2*cr+1)),fr+(1:(2*cr+1)))=fspecial('disk',cr);
C = normxcorr2(h_ext,canvas);
cell_corr=C(1+2*(cr+fr):(size(Img_in,1)+2*(cr+fr)),1+2*(cr+fr):(size(Img_in,2)+2*(cr+fr)));
