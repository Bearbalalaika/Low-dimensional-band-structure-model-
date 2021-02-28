function C =tetra_profile(x, y)
[X, Y] = meshgrid(x, y);
a = 1;

d = sqrt(6)/3*a;

C1x = d/(sqrt(3)*a/3);
C1y = d/(a/3);
C1z = d/(sqrt(6)*a/3);
%P1 = C1x*X + C1y*Y + C1z*Z -d;

C2x = d/(sqrt(3)*a/3);
C2y = -d/(a/3);
C2z = d/(sqrt(6)*a/3);
%P2 = C2x*X + C2y*Y + C2z*Z -d;


C3x = d/(-sqrt(3)*a/6);
C3y = 0;
C3z = d/(sqrt(6)*a/3);

Ca = (C1x*X + C1y*Y)/(C1z);
Cb = (C2x*X + C2y*Y)/(C2z);
Cc = (C3x*X + C3y*Y)/(C3z);

Ca(find((Ca<=Cb)|(Ca<=Cc))) = 0;
Cb(find((Cb<=Ca)|(Cb<=Cc))) = 0;
Cc(find((Cc<=Ca)|(Cc<=Cb))) = 0;



C = Ca+Cb+Cc;