function prct=cheap_prctil(input,perc)
input=reshape(input,[prod(size(input)) 1]);
[B,IX] = sort(input);
if ceil((prod(size(input))+1)*(perc/100))>numel(B)
    prct=B(end);
else
prct=B(ceil((prod(size(input))+1)*(perc/100)));
end