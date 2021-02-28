function M = Dy(X, Y, Z, sym)
    if nargin < 4; sym = 0; end;
    [Sy, Sx, Sz] = size(Y);
    Stot = Sx*Sy*Sz;
    dy = Y(2:Sy, :, :)-Y(1:(Sy-1), :, :);    
    a = 0*Y; b = 0*Y; c = 0*Y;
    
    dyp =  dy;
    dyp(Sy,:,:) = dy(Sy-1,:,:);

    dyn =  dy; i=2:Sy;
    dyn(i,:,:) = dy(i-1,:,:);
    
    a = -dyn./(2*dyp);
    b = dyp./(2*dyn);

    c = -a-b;
        
    alfa = a; beta = b; gamma = c;

    %Boundary condition
    alfa(Sy,:,:) = 0;
    beta(1,:,:) = 0;
                         
    A = reshape(alfa, Stot, 1); B = reshape(beta, Stot, 1); C = reshape(gamma, Stot, 1); 
    M = spdiags([A C B], [-1 0 1], Stot, Stot);
    if sym M = (invLy(X, Y, Z).^2)*M;
    else M = invLy(X, Y, Z)*M*invLy(X, Y, Z); end