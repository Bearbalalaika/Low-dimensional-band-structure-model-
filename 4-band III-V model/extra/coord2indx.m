function ij = coord2indx(x, y, x0, y0)
ij = round([interp1(x,1:length(x),  x0) interp1(y,1:length(y),  y0)]);
