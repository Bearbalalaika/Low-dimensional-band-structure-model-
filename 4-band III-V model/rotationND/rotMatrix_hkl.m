%A = rotMatrix(hkl) where hkl is a vector containing the miller indicies
%for rotated z, y = [-kh0]
%See Review Strain and strain relaxation in semiconductors
%JOURNAL OF MATERIALS SCIENCE: MATERIALS IN ELECTRONICS 8 (1997) 337-375
%
%For the function not used:
%See fig in J. Phys. D: Appl. Phys., 1971, Vol. 4. (J. TURLEY and G. SINES pp 1731) 

function R = rotMatrix(hkl)
    
hkl = hkl./(sqrt(hkl(1)^2 + hkl(2)^2 + hkl(3)^2));
h = hkl(1); k = hkl(2); l = hkl(3);

if (hkl == [0 0 1])
    R = [1 1 0; -1 1 0; 0 0 sqrt(2)]/sqrt(2);
else
    
sq = sqrt(1-l^2);

%In Steps  ORIGINAL y = [-kh0]
%XY-rot
%R1 = [h/sq k/sq 0;
%     -k/sq h/sq 0;
%      0      0  1];
%Z-rot
%R2 = [l 0 -sq;
%      0 1 0;
%      sq 0 l];
%R = R2*R1;

%Full ORIGINAL y = [-kh0]
R = [h*l/sq k*l/sq -sq;
     -k/sq   h/sq    0;
       h      k       l];
end


%THIS FUNCTION SHOULE BE CHECKED BUT IS MORE GENERAL
function R0 = rotM(th, hkl)
%Alternatively use euler angles a and b
h = hkl(1); k = hkl(2); l = hkl(3);
sin_a = k/sqrt(h^2 + k^2);
cos_a = h/sqrt(h^2 + k^2);
sin_b = l/sqrt(h^2 + k^2 + l^2);
cos_b = sqrt((h^2+k^2))/sqrt(h^2 + k^2 + l^2);
  
A = cos_a*cos_b;
B = sin_a*cos_b;
C = sin_b;
D = -cos_a*sin_b;
E = -sin_a;
F = -sin_a*sin_b;
G = cos_a;
H = cos_b;

R0 = [A                               B                    C;
     D*sin(th)+E*cos(th)     F*sin(th)+G*cos(th)         H*sin(th);
     D*cos(th)-E*sin(th)     F*cos(th)-G*sin(th)         H*cos(th);];
L = [0 0 1; 0 1 0; 1 0 0];
R0 = (L)*R0*inv(L);
