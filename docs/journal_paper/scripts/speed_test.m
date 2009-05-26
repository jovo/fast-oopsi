CC=C(:);
CCC=reshape(CC,Nc,Meta.T);
% bb=b(:);
for t=1:100
D   = F-P.a*(CCC+b);
D   = F-P.a*(reshape(CC,Nc,Meta.T)+b);
end