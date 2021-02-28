function M = B_k2y(f, X, Y, Z, Ay)

[Sy, Sx, Sz] = size(X); Stot = Sx*Sy*Sz;
    
F = spdiags(reshape(f, Stot, 1), 0, Stot, Stot);
Ay = spdiags(reshape(Ay, Stot, 1), 0, Stot, Stot);

M = k2y(f, X, Y, Z) + ky(X, Y, Z)*F*Ay + Ay*F*ky(X, Y, Z) + Ay*F*Ay;
