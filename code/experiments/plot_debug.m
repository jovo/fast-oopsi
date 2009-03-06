fnum=fnum+1; figure(fnum), clf,

mn(1) = min([min(P.a) min(P.b)]);
mx(1) = (max([max(P.a) max(P.b)])-mn(1))/60;

for q=qs
    mn(1+q) = min(Phat{q}.a(:));
    mx(1+q) = (max(Phat{q}.a(:))-mn(1))/60;
end

mnn=min(mn);
mxx=max(mn);

nrows=numel(qs);
ncols=3;

for q=qs
    subplot(nrows,ncols,1+(q-1)*ncols),
    Pl.ylab=I{q}.label;
    if q==1, Pl.tit='Sum of Filters'; else Pl.tit=[]; end
    Plot_im(Pl,sum(I{q}.P.a,2))
end


for i=1:Nc, 
    
    for q=qs
        subplot(nrows,ncols,i+1+(q-1)*ncols),
        Pl.ylab=[];
        if q==1,  Pl.tit=['Filter ' num2str(i)]; else Pl.tit=[]; end
        Plot_im(Pl,(I{q}.P.a(:,i)-mn(q+1))/mx(q+1))
    end
end
set(gca,'XTick',[0:5:25],'YTick',[0:5:25])

% print fig
wh=[7 5];   %width and height
set(fnum,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc',[fig_fname 'est_filers'])