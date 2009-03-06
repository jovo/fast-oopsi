function MakeTif(F,Sim)

F=F-min(F(:));
F=round((2^Sim.bits-1)*F/max(F(:)));
if Sim.bits==8, F=uint8(F); else F=uint16(F); end
for i=1:Sim.T
    if i==1, mod='overwrite'; else mod='append'; end
    imwrite(reshape(F(i,:),Sim.w,Sim.h),Sim.tif_name,'tif','Compression','none','WriteMode',mod)
end
