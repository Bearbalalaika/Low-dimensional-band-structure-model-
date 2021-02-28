    function M = invLz(X, Y, Z)
    [Sy, Sx, Sz] = size(Z);
    Stot = Sx*Sy*Sz;
    if Sz>1
        Stot = Sx*Sy*Sz;
        dzp = Z(:,:, 2:Sz)-Z(:, :,1:(Sz-1)); dzp(:, :, Sz) = dzp(:,:, Sz-1);
        dzn(:,:,1) = dzp(:,:,1); dzn(:,:,2:Sz) = dzp(:,:,1:(Sz-1));
        invL = 1./sqrt((dzp+dzn)/2);
        M = spdiags(reshape(invL, Stot, 1), 0, Stot, Stot);
    else
        M = spdiags(ones(Stot, 1), 0, Stot, Stot);
    end