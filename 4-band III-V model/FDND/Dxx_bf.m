function M = Dxx_bf(Q, X, Y, Z)
    [Sy, Sx, Sz] = size(Y);
    Stot = Sx*Sy*Sz;
    q = spdiags(reshape(Q, Stot, 1), 0, Stot, Stot);
    
    M = ( Dx_f(X, Y, 0)*q*Dx_b(X, Y, Z) + Dx_b(X, Y, 0)*q*Dx_f(X, Y, Z))/2;
    