function [ff_corr]=ffcorr(input,radius,showflag);

if nargin<3
    showflag=0;
end
input=double(input);
h = fspecial('disk',radius);
input_filt=imfilter(double(input), h, 'replicate');
ff_corr=real(input./input_filt);

if showflag
figure, imshow(input,[min(min(input)) max(max(input))]), title('before ff corr')   
figure, imshow(ff_corr,[min(min(ff_corr)) max(max(ff_corr))]), title('after ff corr')
end