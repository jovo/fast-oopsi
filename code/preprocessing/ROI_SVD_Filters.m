Phat{q}=P;

for i=1:Nc
    COLS=COLS1{i};
    ROWS=ROWS1{i};
    if i==Nc, COLS(2)=Ncols_ROI; ROWS(2)=Nrows_ROI; end
    ind=[];
    for j=COLS(1):COLS(2)
        ind=[ind (j-1)*Nrows_ROI + [ROWS(1):ROWS(2)]];
    end
    FF          = F(:,ind);
    [U,S,V]     = svd(FF,0);
    Phat{q}.a(:,i)   = 0;
    if sum(V(:,1))<0, V(:,1)=-V(:,1); end
    Phat{q}.a(ind,i) = V(:,1);
end

mn = min([min(a_b) min(P.a(:)) min(Phat{q}.a(:))]);
mx = (max([max(a_b) max(P.a(:)) max(Phat{q}.a(:))])-mn)/60;

fnum=fnum+1; fig=figure(fnum); clf,ncols=1+Nc;
subplot(1,ncols,1), image(reshape((a_b-mn)/mx,Nrows_ROI,Ncols_ROI))%, colorbar
for i=1:Nc
    subplot(1,ncols,1+i), image(reshape((Phat{q}.a(:,i)-mn)/mx,Nrows_ROI,Ncols_ROI))%, colorbar
    title(['SVD' num2str(i)])
end
wh=[7 2];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-deps','Multi_ROI_SVD')
I{q}.label    = 'ROI SVD Filter';
