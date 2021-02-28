
function C_QD = Side_QWR(x, y, z,r,l, alcont,core, asym,gauss)
 
r1=30;
r2=70;
 
 
 [X1, Y1, Z1] = meshgrid(x, y, z);

 %R1 = sqrt(X1.^2/asym + asym*Y1.^2+(Z1).^2);
 
 cd = ones(size(Y1));
 
 cd=cd.*0.4;
 % Plane with conctant al conc
 

 cd(find(abs(Y1)<=(l/2))) = alcont;
 
 
 % GaAs QWR inside the plane
 
% R2 = sqrt(X1.^2/asym + asym*Y1.^2);
 
  R1 = sqrt(  asym*Y1.^2+(X1-(r1/2-r)).^2);

  
 R2 = sqrt( asym*Y1.^2+(X1-(r2/2)).^2);
 
 
 cd(find(((R1<=r1/2)&R2>r2/2))) = core;
 
 
 cd=cd.*exp(-(X1.^2+Y1.^2+Z1.^2)*gauss);
 
 
 
 C_QD=cd;


