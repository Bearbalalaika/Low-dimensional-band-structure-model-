%IB = broaden(E, I, FWHM)
%E = Energy
%I = Intensity
%FWHM = bradening parameters [Gaussian Lorentzian]x = 

function IB = broaden(E, I, FWHM)

if (nargin == 2)
   FWHM = I;
   I = E;
   E = 1:length(I);
end

dE = mean(diff(E));
c = round(mean(E));

if length(FWHM) == 1
    FWHM = [FWHM(1) 0];
end

B = analvoigt(1, c, FWHM(2), FWHM(1), E);


B = B/trapz(B);

lengthI = length(I); I = [ones(1, length(B)+1)*I(1) I ones(1, length(B)+1)*I(length(I))];

IB = conv(I, B);

IB(1:(length(B)+1)) = []; IB = IB(1:(lengthI + length(B)));

IB(1:floor(interp1(E, 1:length(E), c)-dE/2)) = [];
IB = IB(1:length(E));