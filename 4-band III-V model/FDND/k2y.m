function M = k2y(F, X, Y, Z)
    [Sy, Sx, Sz] = size(Y);
    Stot = Sx*Sy*Sz;
    dy = Y(2:Sy, :, :)-Y(1:(Sy-1), :, :);    
    a = 0*Y; b = 0*Y; c = 0*Y;
    
    %Mean mass function F
    Fn = F; Fp = F;
    i = 1:(Sy-1); Fp(i,:,:) = (Fp(i,:,:) + Fp(i+1,:,:))/2; 
    i = 2:(Sy); Fn(i,:,:) = (Fn(i,:,:) + Fn(i-1,:,:))/2;
    
    dyp =  dy;
    dyp(Sy,:,:) = dy(Sy-1,:,:);

    dyn =  dy; i=2:Sy;
    dyn(i,:,:) = dy(i-1,:,:);
    
    a = Fp./dyp;
    b = Fn./dyn;

    c = -a-b;
        
    alfa = a; beta = b; gamma = c;

    %Boundary condition
    alfa(Sy,:,:) = 0;
    beta(1,:,:) = 0;
                         
    A = reshape(alfa, Stot, 1); B = reshape(beta, Stot, 1); C = reshape(gamma, Stot, 1);
    M = spdiags([A C B], [-1 0 1], Stot, Stot);
    M = invLy(X, Y, Z)*M*invLy(X, Y, Z);
    M = (-1i)^2*M;