function M = D2z(F, X, Y, Z)
    [Sy, Sx, Sz] = size(Z);
    Stot = Sx*Sy*Sz;
    dz = Z(:, :, 2:Sz)-Z(:, :, 1:(Sz-1));    
    a = 0*X; b = 0*X; c = 0*X;
    
    %Mean mass function F
    Fn = F; Fp = F;
    i = 1:(Sz-1); Fp(:, :,i) = (Fp(:, :,i) + Fp(:, :,i+1))/2; 
    i = 2:(Sz); Fn(:, :,i) = (Fn(:, :,i) + Fn(:, :,i-1))/2;
    
    dzp = dz;
    dzp(:,:,Sz) = dz(:,:,Sz-1);

    dzn = dz; i=2:Sz;
    dzn(:,:,i) = dz(:,:,i-1);
    
    a = Fp./dzp;
    b = Fn./dzn;

    c = -a-b;
    
    alfa = a; beta = b; gamma = c;

    
    %Boundary condition
    alfa(:,:, Sz) = 0;
    beta(:,:,1) = 0;
                         
    A = reshape(alfa, Stot, 1); B = reshape(beta, Stot, 1); C = reshape(gamma, Stot, 1); 
    M = spdiags([A C B], [-Sx*Sy 0 Sx*Sy], Stot, Stot);
    M = invLz(X, Y, Z)*M*invLz(X, Y, Z);
    M = (-1i)^2*M;
    