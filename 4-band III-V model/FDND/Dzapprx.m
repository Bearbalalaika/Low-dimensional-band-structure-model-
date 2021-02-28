function M = Dzapprx(X, Y, Z)
   [Sy, Sx, Sz] = size(Z);
    Stot = Sx*Sy*Sz;
    
    a = -1/2*ones(size(Z));
    b =1/2*ones(size(Z));
    c = -a-b;
    
    alfa = a; beta = b; gamma = c;
    
    %Boundary condition
    alfa(:,:, Sz) = 0;
    beta(:,:,1) = 0;
                         
    A = reshape(alfa, Stot, 1); B = reshape(beta, Stot, 1); C = reshape(gamma, Stot, 1); 
    M = spdiags([A C B], [-Sx*Sy 0 Sx*Sy], Stot, Stot);
    M = invLz(X, Y, Z)*M*invLz(X, Y, Z);
