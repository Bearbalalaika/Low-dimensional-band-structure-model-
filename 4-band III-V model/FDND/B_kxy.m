function M = B_kxy(f, X, Y, Z, Ax, Ay)

    [Sy, Sx, Sz] = size(X); Stot = Sx*Sy*Sz;
     F = spdiags(reshape(f, Stot, 1), 0, Stot, Stot);
    Ax = spdiags(reshape(Ax, Stot, 1), 0, Stot, Stot);
    Ay = spdiags(reshape(Ay, Stot, 1), 0, Stot, Stot);

    
    M =     kx(X, Y, Z)*F*ky(X, Y, Z) + kx(X, Y, Z)*F*Ay + Ax*F*ky(X, Y, Z) + Ax*F*Ay;
    M = M + ky(X, Y, Z)*F*kx(X, Y, Z) + ky(X, Y, Z)*F*Ax + Ay*F*kx(X, Y, Z) + Ay*F*Ax;
    M = M / 2;
    