%s.x0 = 0; s.y0 = 0; s.z0 = 0; s.r = 20; s.h = 5; C = truncsphere(s);
function C = truncsphere(x, y, z, shape)
if nargin<3
    shape = x;
    x = -20:0.4:20;
    y = -20:0.4:20;
    z = -20:0.4:20;
end
   
    C = buildTSPHERE(x, y, z, shape.x0, shape.y0, shape.z0, shape.r, shape.h);

function C = buildTSPHERE(x, y, z, x0, y0, z0, r, h)
[X, Y, Z] = meshgrid(x-x0, y-y0, z-z0);

Rmax = (r^2+h^2)/(2*h);
R = sqrt(X.^2+Y.^2+(Z+(Rmax-h)).^2);
C = zeros(size(X));

C(find((R<Rmax)&(Z>=0))) = 1;


