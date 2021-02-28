
function C_QD = QD_lens(x, y, z, a1,a2, l, asym,gauss, C_QD)
 

 
 [X1, Y1, Z1] = meshgrid(x, y, z);

 %clear x y z
 R1 = sqrt(X1.^2/asym + asym*Y1.^2+(Z1-(a1/2-l)).^2);
 cd = zeros(size(R1));
 %c(find(   (R<a/2)&(Z<=(l/2))&(Z>(-l/2))       )) = 1;
 R2 = sqrt(X1.^2/asym + asym*Y1.^2+(Z1-(a2/2)).^2);
 
 
 cd(find(((R1<=a1/2)&R2>a2/2))) = 1;
 
 
 cd=cd.*exp(-(X1.^2+Y1.^2+Z1.^2)*gauss);
 
 cd=1-cd;
 
 
 C_QD=C_QD.*cd;