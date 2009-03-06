[U,S,V] = svd(F-repmat(mean(F),Sim.T,1),0);
PCs = 1:3;
Denoised = U(:,PCs)*S(PCs,PCs)*V(:,PCs)';

if MakMov==1
    for i=1:Sim.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(Denoised(i,:),width,height),'Multi_Mov_SVD_no_mean.tif','tif','Compression','none','WriteMode',mod)
    end
end

a_b=a_b/norm(a_b);
mn = min([min(a_b) min(min(V(:,1:Nc+1)))]);
mx = (max([max(a_b) max(max(V(:,1:Nc+1)))])-mn)/60;

fnum=fnum+1; fig=figure(fnum); clf,ncols=2+Nc;
subplot(1,ncols,1), image(reshape((a_b-mn)/mx,width,height))%, colorbar
title('mean subtracted svds')
Phat{q}=P;
for i=1:Nc
    if sum(V(:,1))<0, V(:,1)=-V(:,1); end
    Phat{q}.a(:,i) = V(:,i);
    if max(F*Phat{q}.a(:,i))<0, Phat{q}.a(:,i)=-Phat{q}.a(:,i); end
    subplot(1,ncols,1+i), image(reshape((Phat{q}.a(:,i)-mn)/mx,width,height))%, colorbar
    title(['SVD' num2str(i)])
end
% Phat{q}.b = V(:,Nc+1);
subplot(1,ncols,ncols), image(reshape((V(:,Nc+1)-mn)/mx,width,height))%, colorbar
title(['SVD' num2str(i+1)])
wh=[7 2];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-deps','Multi_SVD_no_mean')
I{q}.label    = 'SVD Filter (no mean)';
