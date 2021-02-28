function M = Dy_f(X, Y, Z)
    [Sy, Sx, Sz] = size(Y);
    Stot = Sx*Sy*Sz;
    dy = Y(2, 1, 1)-Y(1, 1, 1);    
    C = ones(Sy,Sx, Sz);
    N = ones(Sy,Sx, Sz);

    N(1,:,:) = 0; 
    
    C = reshape(C, Stot, 1);
    N = reshape(N, Stot, 1);
    
    M = spdiags([-C N], [0 1], Stot, Stot)/dy;