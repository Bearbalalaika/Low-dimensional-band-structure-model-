function M = Dz(X, Y, Z, sym)
    if nargin < 4; sym = 0; end;
   [Sy, Sx, Sz] = size(Z);
    Stot = Sx*Sy*Sz;
    dz = Z(:, :, 2:Sz)-Z(:, :, 1:(Sz-1));    
    a = 0*X; b = 0*X; c = 0*X;
       
    dzp = dz;
    dzp(:,:,Sz) = dz(:,:,Sz-1);

    dzn = dz; i=2:Sz;
    dzn(:,:,i) = dz(:,:,i-1);
    
    a = -dzn./(2*dzp);
    b = dzp./(2*dzn);

    c = -a-b;
    
    alfa = a; beta = b; gamma = c;

    
    %Boundary condition
    alfa(:,:, Sz) = 0;
    beta(:,:,1) = 0;
                         
    A = reshape(alfa, Stot, 1); B = reshape(beta, Stot, 1); C = reshape(gamma, Stot, 1); 
    M = spdiags([A C B], [-Sx*Sy 0 Sx*Sy], Stot, Stot);
    if sym M = M*invLz(X, Y, Z).^2;
    else M = invLz(X, Y, Z)*M*invLz(X, Y, Z); end
