[U,S,V] = svd(F,0);
PCs = 1:Nc;
Denoised = U(:,PCs)*S(PCs,PCs)*V(:,PCs)';

if MakMov==1
    for i=1:Sim.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(Denoised(i,:),width,height),'Multi_Mov_Denoised.tif','tif','Compression','none','WriteMode',mod)
    end
end

a_b=a_b/norm(a_b);
mn = min([min(a_b) min(min(V(:,1:Nc+1)))]);
mx = (max([max(a_b) max(max(V(:,1:Nc+1)))])-mn)/60;

fnum=fnum+1; fig=figure(fnum); clf,ncols=2+Nc;
subplot(1,ncols,1), image(reshape((a_b-mn)/mx,width,height))%, colorbar
title('svds')
Phat{q}=P;
for i=1:Nc
    if sum(V(:,i))<0, V(:,i)=-V(:,i); end
    Phat{q}.a(:,i) = V(:,i);
    subplot(1,ncols,1+i), image(reshape((Phat{q}.a(:,i)-mn)/mx,width,height))%, colorbar
    title(['SVD' num2str(i)])
end
% Phat{q}.b = V(:,Nc+1);
subplot(1,ncols,ncols), image(reshape((V(:,Nc+1)-mn)/mx,width,height))%, colorbar
title(['SVD' num2str(i+1)])
wh=[7 2];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-deps','Multi_SVD')
I{q}.label    = 'SVD Filter';
