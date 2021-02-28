%s.x0 = 0; s.y0 = 0; s.z0 = 0; s.a = 20; s.b = 20; s.h = 5; s.f = 1; C = truncpyr(s);
function C = truncpyr(x, y, z, shape)
if nargin<3
    shape = x;
    x = -20:0.5:20;
    y = -20:0.5:20;
    z = -20:0.5:20;
end
    c = 0;
    if isfield(shape, 'circularity'); c = shape.circularity; end
    C = buildTPYR(x, y, z, shape.x0, shape.y0, shape.z0, shape.a, shape.b, shape.h, shape.f, c);

function C = buildTPYR(x, y, z, x0, y0, z0, a, b, h, f, c)
[X, Y, Z] = meshgrid(x-x0, y-y0, z-z0);

X0=sqrt(X.^2 + c*Y.^2);
Y0=sqrt(c*X.^2 + Y.^2);


C = zeros(size(X));

Cx = a*(1-Z/h)/2;
Cy = b*(1-Z/h)/2;

C(find( (abs(X0)<=Cx)&(abs(Y0)<=Cy)&(Z<=(h*f))&(Z>=0) )) = 1; disp('EQUAL')
%C(find( (abs(X0)<Cx)&(abs(Y0)<Cy)&(Z<(h*f))&(Z>0) )) = 1;


