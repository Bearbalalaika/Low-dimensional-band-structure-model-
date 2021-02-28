    
function M = Dyx(F, X, Y, Z)
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    F = spdiags(reshape(F, Stot, 1), 0, Stot, Stot);
    M = Dyapprx(X, Y, Z)*F*Dxapprx(X, Y, Z);