b=rand(100,1);
Z=zeros(10,1);
for i=1:1e3
    ass1=(1+Z)*b';
    ass2=repmat(b,1,10)';
end

for i=1:1e3
    ass3=(1+Z)*b';
    ass4=repmat(b,10,1);
end

