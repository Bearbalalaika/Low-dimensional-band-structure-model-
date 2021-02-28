function M = B_k2z(f, X, Y, Z, Az)

[Sy, Sx, Sz] = size(X); Stot = Sx*Sy*Sz;
    
F = spdiags(reshape(f, Stot, 1), 0, Stot, Stot);
Az = spdiags(reshape(Az, Stot, 1), 0, Stot, Stot);

M = k2z(f, X, Y, Z) + kz(X, Y, Z)*F*Az + Az*F*kz(X, Y, Z) + Az*F*Az;
