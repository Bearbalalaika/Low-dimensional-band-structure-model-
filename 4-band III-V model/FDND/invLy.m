function M = invLy(X, Y, Z)
    [Sy, Sx, Sz] = size(Y);
    Stot = Sx*Sy*Sz;
    if Sy>1
        Stot = Sx*Sy*Sz;
        dyp = Y(2:Sy, :, :)-Y(1:(Sy-1), :, :); dyp(Sy, :, :) = dyp(Sy-1,:,:);   
        dyn(1,:, :) = dyp(1,:, :); dyn((2:Sy),:, :) = dyp((1:(Sy-1)), :, :);
        dyn = reshape(dyn, Sy, Sx, Sz); %Varför ????????????????????
        invL = 1./sqrt((dyp+dyn)/2);
        M = spdiags(reshape(invL, Stot, 1), 0, Stot, Stot);
    else
        M = spdiags(ones(Stot, 1), 0, Stot, Stot);
    end