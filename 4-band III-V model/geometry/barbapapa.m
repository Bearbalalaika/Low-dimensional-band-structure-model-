function C = barbapapa(x, y, z)

[X, Y, Z] = meshgrid(x, y, z);
C = zeros(size(X));

for i = 1:14;
    x = 5*cos(i*pi/10);
    y = 5*sin(i*pi/10);
    
    C = C + sphere(X, Y, Z, x, y, -18+2*i, 1+i/1.5);
end

C(find(C>0)) = 1;


function C = sphere(X, Y, Z, x, y, z, r)
  R = sqrt((X-x).^2+(Y-y).^2+(Z-z).^2);
  C = zeros(size(R));
  C(find(R<r)) = 1;

