function M = invLx(X, Y, Z)
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    if Sx>1
        Stot = Sx*Sy*Sz;
        dxp = X(:, 2:Sx, :)-X(:, 1:(Sx-1), :);
        dxp(:,Sx,:) = dxp(:,Sx-1,:);
        dxn(:,2:Sx,:) = dxp(:,1:(Sx-1), :);
        dxn(:,1,:) = dxp(:,1,:);  
        invL = 1./sqrt((dxp+dxn)/2);
        M = spdiags(reshape(invL, Stot, 1), 0, Stot, Stot);
    else
        M = spdiags(ones(Stot, 1), 0, Stot, Stot);
    end