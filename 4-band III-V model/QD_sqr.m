function C_QD = QD_sqr(x, y, z, asym,gauss, C_QD, alcont)
 

 
 [X1, Y1, Z1] = meshgrid(x, y, z);

 %clear x y z
 R1 = sqrt(X1.^2/asym + asym*Y1.^2+(Z1-(a1/2-l)).^2);
 %cd = zeros(size(R1));
 %c(find(   (R<a/2)&(Z<=(l/2))&(Z>(-l/2))       )) = 1;
 
 
 
 C_QD(find(((z<=20)&(x<=20)&(y<=20)))) = alcont;
 
 
% cd=cd.*exp(-(X1.^2+Y1.^2+Z1.^2)*gauss);
 
 %cd=1-cd;
 
 
 %C_QD=C_QD.*cd;