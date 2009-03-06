function b=ownmedfilt2tst(Img_in,filtsz);

a=Img_in;
[xsz ysz]=size(Img_in);
b=a;
h = waitbar(0,['Applying median filter, radius ' num2str(filtsz)]); 
for i=1:xsz
    if ~mod(i,filtsz)
    waitbar(((i-1)*ysz)/(xsz*ysz),h)
    end
    for ii=1:ysz
        
        %center=(ii-1)*xsz+i;
        %r=zeros(2*filtsz+1);
        
        if (i-filtsz)<1
            xrange=(1:(i-1)+filtsz)';
        elseif (i+filtsz)>xsz
            xrange=(i-filtsz:xsz)';
        else
            xrange=(i-filtsz:(i-1)+filtsz)';
        end
        
         if (ii-filtsz)<1
            yrange=(1:(ii-1)+filtsz)';
        elseif (ii+filtsz)>ysz
            yrange=(ii-filtsz:ysz)';
        else
            yrange=(ii-filtsz:(ii-1)+filtsz)';
         end
        
        yrangemat=yrange(:,ones(length(xrange),1));
        xrangemat=xrange(:,ones(length(yrange),1));
        %r=zeros(length(xrange),length(yrange));
        %{
        if (i-filtsz)>0 & (i+filtsz)<=xsz & ...
           (ii-filtsz)>0 & (ii+filtsz)<=ysz     
                
        for j=(center-filtsz):(center+filtsz)
        r(j-(center-filtsz)+1,:)=j+((-filtsz:filtsz)*xsz);
        end
        range=r(r>0&r<=(xsz*ysz));
        %}
           
        tmp=a(xrangemat,yrangemat');
        b(i,ii)=median(tmp(:));
        end
    end
close(h)
