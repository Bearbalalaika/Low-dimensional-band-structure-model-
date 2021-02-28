function [L1, L2, L3] = luttingerParameters(conc, mparam, LINorINV)
    x_In = imag(conc); x_Al = real(conc);
    if strcmp(LINorINV, 'INV')
        display('Inverse interpolation of Luttinger parameters')
        L1 = 1./((1-x_In-x_Al)/mparam.lpGa(1) + x_In/mparam.lpIn(1) + x_Al/mparam.lpAl(1));
        L2 = 1./((1-x_In-x_Al)/mparam.lpGa(2) + x_In/mparam.lpIn(2) + x_Al/mparam.lpAl(2));
        L3 = 1./((1-x_In-x_Al)/mparam.lpGa(3) + x_In/mparam.lpIn(3) + x_Al/mparam.lpAl(3));
    elseif strcmp(LINorINV, 'LIN')
        display('Linear interpolation of Luttinger parameters')
        L1 = (1-x_In-x_Al)*mparam.lpGa(1) + x_In*mparam.lpIn(1) + x_Al*mparam.lpAl(1);
        L2 = (1-x_In-x_Al)*mparam.lpGa(2) + x_In*mparam.lpIn(2) + x_Al*mparam.lpAl(2);
        L3 = (1-x_In-x_Al)*mparam.lpGa(3) + x_In*mparam.lpIn(3) + x_Al*mparam.lpAl(3);
    else
        L1 = NaN; L2 = NaN; L3 = NaN;
    end