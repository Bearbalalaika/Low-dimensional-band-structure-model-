function plotDirection(hkl, type)
if nargin<2
    type = 'k';
end

hold on;
%R = rotMatrix(0, hkl);
%x0 = [1 0 0]'; x = R*x0
hkl = 15*hkl/sqrt(hkl(1)^2 + hkl(2)^2 + hkl(3)^2);
x = hkl;
h=plot3([0 x(1)], [0 x(2)], [0 x(3)], type); set(h, 'linewidth', 3);

axis(15*[-1 1 -1 1 -1 1]); daspect([1 1 1]); box on;