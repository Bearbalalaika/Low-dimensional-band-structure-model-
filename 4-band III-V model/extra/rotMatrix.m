function out = rotMatrix(M, angle, imethod ,xc, yc)
SZE = size(M);
out = zeros(SZE);

if nargin < 3
   imethod = 'linear';
end

if nargin < 4
   xc = SZE(2)/2 + 0.5;
end

if nargin <5
   yc = SZE(1)/2 + 0.5;
end


i = 1:SZE(2);
for   j = 1:SZE(1);

X = sin(angle) * (j - yc) + cos(angle)* (i-xc);
Y = cos(angle)*(j - yc) -  sin(angle) * (i-xc);

XI(j,:) = X;
YI(j,:) = Y;
end
out(i,:) = interp2(M, xc+XI, yc+YI, imethod);
