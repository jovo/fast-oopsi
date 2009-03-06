function L_out=breakcells(L_in);

L_out=L_in.*logical(watershed(imcomplement(bwdist(~L_in))));