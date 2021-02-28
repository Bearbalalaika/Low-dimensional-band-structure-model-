function M = kx(X, Y, Z)
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    dx = X(:, 2:Sx, :)-X(:, 1:(Sx-1), :);    
    a = 0*X; b = 0*X; c = 0*X;
    
    a = -1/2*ones(size(Y));
    b = 1/2*ones(size(Y));
    c = -a-b;
        
    alfa = a; beta = b; gamma = c;
    
    %Boundary condition
    alfa(:,Sx,:) = 0;
    beta(:,1,:) = 0;
                         
    A = reshape(alfa, Stot, 1); B = reshape(beta, Stot, 1); C = reshape(gamma, Stot, 1);    
    M = spdiags([A C B], [-Sy 0 Sy], Stot, Stot);
    M = invLx(X, Y, Z)*M*invLx(X, Y, Z);
    M = (-1i)*M;