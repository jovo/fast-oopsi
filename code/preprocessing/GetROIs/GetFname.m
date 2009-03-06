function fname = GetFname(folder, head, t, foot)

if t<10
    pad='000';
elseif t>9 & t<100
    pad='00';
elseif t>99 & t<1000
    pad='0';
else
    pad='';
end
fname=[folder head pad num2str(t) foot];