function M = Dx(X, Y, Z, sym)
    if nargin < 4; sym = 0; end;
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    dx = X(:, 2:Sx, :)-X(:, 1:(Sx-1), :);    
    a = 0*X; b = 0*X; c = 0*X;
    
    dxp =  dx;
    dxp(:,Sx,:) = dx(:,Sx-1,:);

    dxn =  dx; i=2:Sx;
    dxn(:,i,:) = dx(:,i-1,:);
    
    a = -dxn./(2*dxp);
    b = dxp./(2*dxn);

    c = -a-b;
        
    alfa = a; beta = b; gamma = c;
    
    %Boundary condition
    alfa(:,Sx,:) = 0;
    beta(:,1,:) = 0;
                         
    A = reshape(alfa, Stot, 1); B = reshape(beta, Stot, 1); C = reshape(gamma, Stot, 1);    
    M = spdiags([A C B], [-Sy 0 Sy], Stot, Stot);
    if sym M = (invLx(X, Y, Z).^2)*M;
    else M = invLx(X, Y, Z)*M*invLx(X, Y, Z); end
