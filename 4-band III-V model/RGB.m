function C = RGB(hh)

R = 2*hh; i=find(R>1); R(i) = 1; 

B = 2-2*hh; i=find(B>1); B(i) = 1; 

G = (R + B)-1;

C = [R G B];