X       = [C 1+Z];
Phat{q} = P;
for i=1:Npixs
    Y   = F(:,i);
    B   = X\Y;
    for j=1:Nc
        Phat{q}.a(i,j) = B(j);
    end
    Phat{q}.b(i) = B(end);
end
fnum=1; fig=figure(fnum); clf, ncols=2+Nc;
mn = min([min(a_b) min(Phat{q}.a(:)) min(Phat{q}.b(:))]);
mx = (max([max(a_b) max(Phat{q}.a(:)) max(Phat{q}.b(:))])-mn)/60;

subplot(1,ncols,1), image(reshape((a_b-mn)/mx,Nrows_ROI,Ncols_ROI))%, colorbar
for i=1:Nc
    subplot(1,ncols,1+i), image(reshape((Phat{q}.a(:,i)-mn)/mx,Nrows_ROI,Ncols_ROI))%, colorbar
    title(['a' num2str(i)])
end
subplot(1,ncols,ncols), image(reshape((Phat{q}.b-mn)/mx,Nrows_ROI,Ncols_ROI)), title('b')%, colorbar
wh=[7 2];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-deps','Multi_Filters')
I{q}.label    = 'Spike Filter';
