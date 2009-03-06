function [L_out]=regionpropscheck(L,ar_lowth,ar_highth,ecc,ext,smooth);

if ~exist('ecc') || prod(size(ecc))==0; ecc=1;end
if ~exist('ext') || prod(size(ext))==0; ext=0;end
if ~exist('ar_lowth') || prod(size(ar_lowth))==0; ar_lowth=0;end
if ~exist('ar_highth') || prod(size(ar_highth))==0; ar_highth=Inf;end
if exist('smooth') && prod(size(smooth))~=0; 
    %L=medfilt2(L,[2*smooth+1 2*smooth+1]);
    L=ownmedfilt2(L,smooth);
end
props=regionprops(bwlabel(L),'Area','Eccentricity','Extent');
idx=find([props.Area]<=ar_highth&[props.Area]>=ar_lowth...
    &[props.Eccentricity]<=ecc&[props.Extent]>=ext);
L_out=ismember(bwlabel(L),idx);


