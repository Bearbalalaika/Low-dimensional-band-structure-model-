function C = buildTriangle(x, y, z, a, x0, y0)
    [X, Y, Z] = meshgrid(x-x0, y-y0, z);
    L = a;
    C = zeros(size(X));
    F = a/2-1/sqrt(3)*(X+a/4);
    C(find(abs(Y)<abs(F))) = 1;
    C(find(X<(-a/4))) = 0;
    C(find(X>( (sqrt(3)/2-1/4))*a  )) = 0;
    