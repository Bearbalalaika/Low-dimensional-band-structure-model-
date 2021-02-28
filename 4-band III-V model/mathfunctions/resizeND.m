function  A = resizeND(A, sx, sy, sz)

[Sy, Sx, Sz] = size(A);

method = 'linear';
%method = 'cubic';

Gx = 0; Gy = 0; Gz = 0;

if (nargin==2)
    if ((Sx/sx)>1) Gx = 0.85*Sx/sx; end
    A = gaussBlur(A, Gx);
    x0 = 1:Sx;
    x = linspace(1, Sx, sx);
    A = interp1(x0, A, x, method);
elseif (nargin==3)
    if ((Sx/sx)>1) Gx = 0.85*Sx/sx; end; if ((Sy/sy)>1) Gy = 0.85*Sy/sy; end
    A = gaussBlur(A, Gx, Gy);
    x0 = 1:Sx; y0 = 1:Sy;
    [X0, Y0] = meshgrid(x0, y0);
    x = linspace(1, Sx, sx);
    y = linspace(1, Sy, sy);
    [X, Y] = meshgrid(x, y);
    A = interp2(X0, Y0, A, X, Y, method);
elseif (nargin==4)
    if ((Sx/sx)>1) Gx = 0.85*Sx/sx; end;
    if ((Sy/sy)>1) Gy = 0.85*Sy/sy; end
    if ((Sz/sz)>1) Gz = 0.85*Sz/sz; end
    A = gaussBlur(A, Gx, Gy, Gz);
    x0 = 1:Sx; y0 = 1:Sy; z0 = 1:Sz;
    [X0, Y0, Z0] = meshgrid(x0, y0, z0);
    x = linspace(1, Sx, sx);
    y = linspace(1, Sy, sy);
    z = linspace(1, Sz, sz);
    [X, Y, Z] = meshgrid(x, y, z);
    A = interp3(X0, Y0, Z0, A, X, Y, Z, method);
end
