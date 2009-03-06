function  h= getAxes(fh,which)
    ha = findall(fh,'Type','Axes');
    ph= get(ha,'Position');
    if ph{1}(1) == 0 
        t{1} = ha(1); t{2} = ha(2);
    else
        t{1} = ha(2); t{2} = ha(1);
    end
    h = t{which};
end