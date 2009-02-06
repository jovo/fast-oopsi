Phat{q}=P;
for i=1:Nc
    COLS=COLS1{i};
    ROWS=ROWS1{i};
    if i==Nc, COLS(2)=Ncols_ROI; ROWS(2)=Nrows_ROI; end
    ind=[];
    for j=COLS(1):COLS(2)
        ind=[ind (j-1)*Nrows_ROI + [ROWS(1):ROWS(2)]];
    end
    FF=0*F;
    FF(:,ind)   = F(:,ind);
    V           = mean(FF);
    if sum(V)<0, V=-V; end
    Phat{q}.a(:,i) = V/norm(V);
end

mn = min([min(a_b) min(P.a(:)) min(Phat{q}.a(:))]);
mx = (max([max(a_b) max(P.a(:)) max(Phat{q}.a(:))])-mn)/60;

fnum=fnum+1; fig=figure(fnum); clf,ncols=1+Nc;
subplot(1,ncols,1), image(reshape((a_b-mn)/mx,Nrows_ROI,Ncols_ROI))%, colorbar
for i=1:Nc
    subplot(1,ncols,1+i), image(reshape((Phat{q}.a(:,i)-mn)/mx,Nrows_ROI,Ncols_ROI))%, colorbar
    title(['mean' num2str(i)])
end
wh=[7 2];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-deps','Multi_ROI_mean')
I{q}.label    = 'ROI mean Filter';
