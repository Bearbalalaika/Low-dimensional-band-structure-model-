function M = B_k2x(f, X, Y, Z, Ax)

[Sy, Sx, Sz] = size(X); Stot = Sx*Sy*Sz;
    
F = spdiags(reshape(f, Stot, 1), 0, Stot, Stot);
Ax = spdiags(reshape(Ax, Stot, 1), 0, Stot, Stot);

M = k2x(f, X, Y, Z) + kx(X, Y, Z)*F*Ax + Ax*F*kx(X, Y, Z) + Ax*F*Ax;
