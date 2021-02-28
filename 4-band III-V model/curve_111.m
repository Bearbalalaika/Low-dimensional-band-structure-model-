%NanoMod 2007-2013
%Author:
%Fredrik Karlsson, PhD, Docent
%Semiconductor Materials, IFM
%Link?ping University, Sweden
%freka@ifm.liu.se
%
%CODE NOT FOR DISTRIBUTION
%COPYRIGHT Fredrik Karlsson

function QD = curve_111(QD, v)
x = 1:size(QD,2); x = x - mean(x);
y = 1:size(QD,1); y = y - mean(y);
z = 1:size(QD,3); z = z - mean(z);

[X, Y, Z] = meshgrid(x, y, z);

    XX = (X-Y)/sqrt(2);
    YY = (X+Y-2*Z)/sqrt(6);
    ZZ = (X+Y+Z)/sqrt(3);
    
r = sqrt(XX.^2 + YY.^2);

ZZ = ZZ-v*r.^2;

R = [1/sqrt(2) -1/sqrt(2)  0;
     1/sqrt(6) +1/sqrt(6) -2/sqrt(6);
     1/sqrt(3) +1/sqrt(3) +1/sqrt(3)];
iR = inv(R);

X2 = iR(1,1)*XX + iR(1,2)*YY + iR(1,3)*ZZ;
Y2 = iR(2,1)*XX + iR(2,2)*YY + iR(2,3)*ZZ;
Z2 = iR(3,1)*XX + iR(3,2)*YY + iR(3,3)*ZZ;


QD = interp3(X, Y, Z, QD, X2, Y2, Z2, 'nearest');
QD(isnan(QD)) = 0;

end
