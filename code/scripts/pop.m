load('../../data/new_foopsi.mat')
%%
fig=figure(21); clf

subplot(131)
hold on;

imagesc(data.image); colormap gray;
for ii =1:length(data.contours),
    patch(data.contours{ii}(:,1),data.contours{ii}(:,2),[1 1 1],'FaceColor','none','EdgeColor',0.3*[1 1 1]);
end;
xlim([0 size(data.image,2)]);
ylim([0 size(data.image,1)]);
%%
cells=[22 29 116 110 16 109]; %16 109 21 76
colors=get(gca,'ColorOrder');

k=0;
for ii=cells
    k=k+1;
    patch(data.contours{ii}(:,1),data.contours{ii}(:,2),[1 1 1],'FaceColor','none','EdgeColor',colors(k,:));    
end

% xlabel('pixel')
% ylabel('pixel')
title('segmented image')
axis('tight')
%%
% [foo idx] = sort(sum(data.fast_oopsi,2));
% siz=size(data.traces);
% for j=1:siz(1);
%     F=data.traces(j,:);
%     F=F-mean(F);
%     F=F/std(F);
%     [h pvals(j)]=kstest(F);
% end
% [sorted_pvals idx] = sort(pvals);

%%
T=600;
xticks=[0:150:T];
xticklabels=round(xticks/30);
yticks=[0.5:1:length(cells)];
yticklabels=ceil(flipud(yticks')');

subplot(132); cla, hold on;
for j=0:length(cells)-1
    plot(z1(data.traces((cells(j+1)),1:T))+j,'Color',colors(j+1,:))
    set(gca,'XTick',xticks,'XTickLabel',[])
end
set(gca,'XTick',xticks,'XTickLabel',xticklabels)
set(gca,'YTick',yticks,'YTickLabel',yticklabels)
axis('tight')
xlabel('time (sec)')
ylabel('cell #')
title('fluorescence')

subplot(133); cla, hold on;
for j=0:length(cells)-1
    plot(z1(data.fast_oopsi((cells(j+1)),1:T))+j,'Color',colors(j+1,:))
    set(gca,'XTick',xticks,'XTickLabel',[])
end
set(gca,'XTick',xticks,'XTickLabel',xticklabels)
set(gca,'YTick',yticks,'YTickLabel',yticklabels)
axis('tight')
xlabel('time (sec)')
% ylabel('cell #')
title('fast filter')
%%
% print fig
name='../../figs/pop';
wh=[7 2];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-dpdf',name)
print('-depsc',name)
saveas(fig,name)