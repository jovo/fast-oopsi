load('../../data/s1m4_20x.mat')

fig=figure(21); 

subplot(131)
hold on; 

imagesc(data.image); colormap gray; 
for ii =1:length(data.contours),
patch(data.contours{ii}(:,1),data.contours{ii}(:,2),[1 1 1],'FaceColor','none','EdgeColor',[1 1 1]); 
end; 
xlim([0 size(data.image,2)]); 
ylim([0 size(data.image,1)]);
xlabel('pixel')
ylabel('pixel')
title('segmented image')

[foo idx] = sort(sum(data.fast_oopsi,2));
% idx=flipud(idx);
T=500;
N=20;
subplot(132)
imagesc(data.traces(idx(end-N:end),1:T)), colormap('gray')
xticks=[0:200:T];
set(gca,'XTick',xticks,'XTickLabel',xticks/10)
xlabel('time (sec)')
ylabel('cell #')
title('fluorescence')

subplot(133)
imagesc(data.fast_oopsi(idx(end-N:end),1:T)), colormap('gray')
set(gca,'XTick',xticks,'XTickLabel',xticks/10)
xlabel('time (sec)')
ylabel('cell #')
title('fast filter')

%% print fig
name='../../figs/pop';
wh=[7 2];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-depsc',name)
print('-dpdf',name)
saveas(fig,name)