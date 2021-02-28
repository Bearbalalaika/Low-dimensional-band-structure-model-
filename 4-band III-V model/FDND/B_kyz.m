function M = B_kyz(f, X, Y, Z, Ay, Az)

    [Sy, Sx, Sz] = size(X); Stot = Sx*Sy*Sz;
     F = spdiags(reshape(f, Stot, 1), 0, Stot, Stot);
    Ay = spdiags(reshape(Ay, Stot, 1), 0, Stot, Stot);
    Az = spdiags(reshape(Az, Stot, 1), 0, Stot, Stot);

    
    M =     ky(X, Y, Z)*F*kz(X, Y, Z) + ky(X, Y, Z)*F*Az + Ay*F*kz(X, Y, Z) + Ay*F*Az;
    M = M + kz(X, Y, Z)*F*ky(X, Y, Z) + kz(X, Y, Z)*F*Ay + Az*F*ky(X, Y, Z) + Az*F*Ay;
    M = M / 2;
    