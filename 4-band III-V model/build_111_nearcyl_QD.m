%NanoMod 2007-2013
%Author:
%Fredrik Karlsson, PhD, Docent
%Semiconductor Materials, IFM
%Link?ping University, Sweden
%freka@ifm.liu.se
%
%CODE NOT FOR DISTRIBUTION
%COPYRIGHT Fredrik Karlsson

function [c] = build_111_nearcyl_QD(x, y, z, a, l, asym,gauss)
 if nargin<6
     asym = 1;
 end
 
 [X0, Y0, Z0] = meshgrid(x, y, z);
 X = (X0-Y0)/sqrt(2);
 Y = (X0+Y0-2*Z0)/sqrt(6);
 Z = (X0+Y0+Z0)/sqrt(3);
 clear X0 Y0 Z0
 R = sqrt(X.^2/asym + asym*Y.^2);
 c = zeros(size(R));
 c(find(   (R<a/2)&(Z<=(l/2))&(Z>(-l/2))       )) = 1;
 gauss
 c=c.*exp(-(X.^2+Y.^2+Z.^2)*gauss);