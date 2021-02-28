function M = Dx_f(X, Y, Z)
    [Sy, Sx, Sz] = size(Y);
    Stot = Sx*Sy*Sz;
    dx = X(1, 2, 1)-X(1, 1, 1);    
    C = ones(Sy,Sx, Sz);
    N = ones(Sy,Sx, Sz);

    N(:,1,:) = 0; 
    
    C = reshape(C, Stot, 1);
    N = reshape(N, Stot, 1);
    
    M = spdiags([-C N], [0 Sy], Stot, Stot)/dx;