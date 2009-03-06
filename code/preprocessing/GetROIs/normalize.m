function y = normalize(x, mini, maxi)

y=mini+(maxi-mini)*(x-min(min(x)))/(max(max(x))-min(min(x)));