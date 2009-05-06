% [x,y,z] = meshgrid(-2:.2:2,-2:.25:2,-2:.16:2);
% v = x.*exp(-x.^2-y.^2-z.^2);
% xslice = [-1.2,.8,2]; yslice = []; zslice = [];
% slice(x,y,z,v,xslice,yslice,zslice)
% colormap hsv

% %%
% figure(9), clf
% for i = -2:.5:2
%     hsp = surf(reshape(F(400,:),width,height));
% %  hsp = surf(linspace(-2,2,20),linspace(-2,2,20),zeros(20)+i);
%  rotate(hsp,[0,1,.8],90)
%  xd = get(hsp,'XData');
%  yd = get(hsp,'YData');
%  zd = get(hsp,'ZData');
%  delete(hsp)
%  slice(x,y,z,v,[-2,2],2,-2) % Draw some volume boundaries
%  hold on
%  slice(x,y,z,v,xd,yd,zd)
% %  hold off
%  axis tight
%  view(-5,10)
%  drawnow
% end

%%
clc
figure(101), clf
subplot(3,1,1)
v=reshape(F',width,height,Sim.T);
vv=v(:,:,100:20:Sim.T);
[XX,YY,ZZ] = meshgrid(x,y,1:size(vv,3));
inc=5;
xslice = []; yslice = []; zslice = 1:inc:inc^2; %100:200:Sim.T;
hh = slice(XX,YY,ZZ,vv,xslice,yslice,zslice);
rotate(hh,[1 0 1],90)
colormap gray
set(gca,'XTick',[],'FontSize',Pl.fs)
set(gca,'YTick',[],'FontSize',Pl.fs)
set(gca,'ZTick',[],'FontSize',Pl.fs)
axis tight
