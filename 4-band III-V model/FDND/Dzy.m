function M = Dzy(F, X, Y, Z)
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    F = spdiags(reshape(F, Stot, 1), 0, Stot, Stot);
    M = Dzapprx(X, Y, Z)*F*Dyapprx(X, Y, Z);

