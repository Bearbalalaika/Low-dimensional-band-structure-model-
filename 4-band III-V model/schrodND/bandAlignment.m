function B = bandAlignment(conc, mparam, VBorCB)
    x_In = imag(conc); x_Al = real(conc);
   
    dEg = polyval(mparam.EgInGaAs, x_In) + polyval(mparam.EgAlGaAs, x_Al)-polyval(mparam.EgInGaAs, 0);
    %dEg = polyval(mparam.EgInGaAs, x_In) + polyval([0.22 1.36 0], x_Al)-polyval(mparam.EgInGaAs, 0);
    if isfield(mparam, 'VBOInGaAs')
        disp(['Determine spatially dependent VB-offset using VBO (see JAP89)']);
        dVB = -(polyval(mparam.VBOInGaAs, x_In)+polyval(mparam.VBOAlGaAs, x_Al)-polyval(mparam.VBOInGaAs, 0)); 
        %dVB=-0.4*dEg;
    elseif isfield(mparam, 'Q')
        disp(['Warning: Global band offset ' num2str(mparam.Q(1))])
        dVB = (1-mparam.Q(1))*dEg;
        %dVB=-0.4*dEg;
        
    end
    dCB = dEg - dVB;
    %dCB = 0.6*dEg;
    %dVB=0.4*dEg;
    if strcmp(VBorCB, 'VB')
        display('Concentration dependent VALENCE BAND potential determined');
        B = dVB;
    elseif strcmp(VBorCB, 'CB')
        display('Concentration dependent CONDUCTION BAND potential determined')
       
       B = dCB;
    else
        B = NaN;
    end