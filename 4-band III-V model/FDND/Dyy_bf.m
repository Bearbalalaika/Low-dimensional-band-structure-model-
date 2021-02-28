function M = Dyy_bf(Q, X, Y, Z)
    [Sy, Sx, Sz] = size(Y);
    Stot = Sx*Sy*Sz;
    q = spdiags(reshape(Q, Stot, 1), 0, Stot, Stot);
    
    M = ( Dy_f(X, Y, 0)*q*Dy_b(X, Y, Z) + Dy_b(X, Y, 0)*q*Dy_f(X, Y, Z))/2;
    