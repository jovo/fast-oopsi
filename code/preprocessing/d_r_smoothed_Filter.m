Phat{q}=P;
Tim.Nc=1;
Phat{q}.a=1;
Phat{q}.b=0;
Phat{q}.gam=P.gam(1);
Phat{q}.lam=P.lam(1);
mn=min(F);
mx=max(F)-mn;
GG=(F-(ones(Sim.T,1)*mn))./(ones(Sim.T,1)*mx);
deconv=0*FF;
reconv=deconv;
for i=1:Npixs
    deconv(:,i) = FOOPSI2_59(GG(:,i),Phat{q},Tim);
    reconv(:,i) = filter(1,[1, -P.gam(1)],deconv(:,i));
end
GG=reconv;
Phat{q}=P;
Tim=Sim;
fnum=fnum+1; fig=figure(fnum); imagesc(reconv')
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-deps','d-r_2D')

fnum=fnum+1; fig=figure(fnum); plot(reconv(:,513))
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-deps','d-r_1D')

if MakMov==1
    for i=1:Sim.T
        if i==1, mod='overwrite'; else mod='append'; end
        imwrite(reshape(reconv(i,:),Nrows_ROI,Ncols_ROI),'Multi_Mov_d-r.tif','tif','Compression','none','WriteMode',mod)
    end
end

I{q}.label='d-r preprocess, then use true filter';
