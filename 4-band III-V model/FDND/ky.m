function M = ky(X, Y, Z)
    [Sy, Sx, Sz] = size(Y);
    Stot = Sx*Sy*Sz;
    dy = Y(2:Sy, :, :)-Y(1:(Sy-1), :, :);    
    a = 0*Y; b = 0*Y; c = 0*Y;
    
    a = -1/2*ones(size(Y));
    b = 1/2*ones(size(Y));
    c = -a-b;
        
    alfa = a; beta = b; gamma = c;

    %Boundary condition
    alfa(Sy,:,:) = 0;
    beta(1,:,:) = 0;
                         
    A = reshape(alfa, Stot, 1); B = reshape(beta, Stot, 1); C = reshape(gamma, Stot, 1); 
    M = spdiags([A C B], [-1 0 1], Stot, Stot);
    M = invLy(X, Y, Z)*M*invLy(X, Y, Z);
    M = (-1i)*M;
    