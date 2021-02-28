function M = Dxy_bf(Q, X, Y, Z)
    [Sy, Sx, Sz] = size(Y);
    Stot = Sx*Sy*Sz;
    q = spdiags(reshape(Q, Stot, 1), 0, Stot, Stot);
    
    M = (Dx_f(X, Y, 0)*q*Dy_f(X, Y, Z) + Dy_f(X, Y, Z)*q*Dx_f(X, Y, Z) + Dx_b(X, Y, Z)*q*Dy_b(X, Y, Z) + Dy_b(X, Y, Z)*q*Dx_b(X, Y, Z))/4;
    