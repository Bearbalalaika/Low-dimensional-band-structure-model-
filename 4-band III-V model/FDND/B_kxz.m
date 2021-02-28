function M = B_kxz(f, X, Y, Z, Ax, Az)

    [Sy, Sx, Sz] = size(X); Stot = Sx*Sy*Sz;
     F = spdiags(reshape(f, Stot, 1), 0, Stot, Stot);
    Ax = spdiags(reshape(Ax, Stot, 1), 0, Stot, Stot);
    Az = spdiags(reshape(Az, Stot, 1), 0, Stot, Stot);

    
    M =     kx(X, Y, Z)*F*kz(X, Y, Z) + kx(X, Y, Z)*F*Az + Ax*F*kz(X, Y, Z) + Ax*F*Az;
    M = M + kz(X, Y, Z)*F*kx(X, Y, Z) + kz(X, Y, Z)*F*Ax + Az*F*kx(X, Y, Z) + Az*F*Ax;
    M = M / 2;
    