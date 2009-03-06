function [bg_subtract, Iobrcbr]=subtractBG(input,radius,showflag);
if nargin<3
    showflag=0;
end
se = strel('disk',radius);
Ie = imerode(input, se);
Iobr = imreconstruct(Ie, input);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
bg_subtract = input-Iobrcbr;
bg_subtract(bg_subtract<0)=0;

if showflag
figure, imshow(input,[min(min(input)) max(max(input))]), title('before BGsubtract')
figure, imshow(Iobrcbr,[min(min(Iobrcbr)) max(max(Iobrcbr))]), title(['Background, radius ' num2str(radius)])
figure, imshow(bg_subtract,[min(min(bg_subtract)) max(max(bg_subtract))]), title('after BGsubtract')
end