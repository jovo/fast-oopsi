function flag = GetFlag(center,deg)
try 
    [x y] = find(center == 1);
    
    if deg == 1
        X = [ones(length(x),1) x];
        Y = [ones(length(y),1) y];
    elseif deg==2
        X = [ones(length(x),1) x x.*x];
        Y = [ones(length(y),1) y y.*y];
    else
        %flag/score the current object
        X = [ones(length(x),1) x x.*x x.*x.*x];
        Y = [ones(length(y),1) y y.*y y.*y.*y];
    end
    [yb,bint,yr,yrint,ystats] = regress(y,X);
    [xb,xbint,xr,xrint,xstats] = regress(x,Y);
    
    flag = min(xstats(4),ystats(4));
catch
    flag = 0.1;
end