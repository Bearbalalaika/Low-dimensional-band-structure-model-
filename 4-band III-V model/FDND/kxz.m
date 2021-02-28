function M = kxz(F, X, Y, Z)
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    F = spdiags(reshape(F, Stot, 1), 0, Stot, Stot);
    M = (Dxapprx(X, Y, Z)*F*Dzapprx(X, Y, Z) + Dzapprx(X, Y, Z)*F*Dxapprx(X, Y, Z))/2;
    M = (-1i)^2*M;
