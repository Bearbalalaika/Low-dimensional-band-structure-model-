function C = rotXYZ_hkl(C, hkl, NORMorINV, bulk, x0, y0, z0)
[sy, sx, sz] = size(C);
if nargin<3
    NORMorINV = 'NORM';
end
if nargin<4
    bulk = C(1, 1, 1);
end
if nargin<5
    x0 = linspace(-1, 1, sx);
end
if nargin<6
    y0 = linspace(-1, 1, sy);
end
if nargin<7
    z0 = linspace(-1, 1, sz);
end


[X0, Y0, Z0] = meshgrid(x0, y0, z0);

Q = rotMatrix_hkl(hkl);

if strcmp(NORMorINV, 'INV') Q = inv(Q); end

X = Q(1, 1)*X0 + Q(1, 2)*Y0 + Q(1,3)*Z0;
Y = Q(2, 1)*X0 + Q(2, 2)*Y0 + Q(2,3)*Z0;
Z = Q(3, 1)*X0 + Q(3, 2)*Y0 + Q(3,3)*Z0;


C = interp3(X0, Y0, Z0, C, X, Y, Z);
C(isnan(C)) = bulk;

