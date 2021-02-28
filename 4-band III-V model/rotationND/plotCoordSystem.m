%plotCoordSystem(rotMatrix(0, [1 1 1]))
%X BLACK  Y RED  Z BLUE
function plotCoordSystem(R);
hold on; 
%%X [001]
x0 = [0 0 1]'; x = R*x0;
h=plot3([0 x0(3)], [0 x0(2)], [0 x0(1)], 'k:'); set(h, 'linewidth', 3);
h=plot3([0 x(3)], [0 x(2)], [0 x(1)], 'k'); set(h, 'linewidth', 3);


%%Y [010]
x0 = [0 1 0]'; x = R*x0;
h=plot3([0 x0(3)], [0 x0(2)], [0 x0(1)], 'r:'); set(h, 'linewidth', 3);
h=plot3([0 x(3)], [0 x(2)], [0 x(1)], 'r'); set(h, 'linewidth', 3);

%%Z [100]
x0 = [1 0 0]'; x = R*x0;
h=plot3([0 x0(3)], [0 x0(2)], [0 x0(1)], 'b:'); set(h, 'linewidth', 3);
h=plot3([0 x(3)], [0 x(2)], [0 x(1)], 'b'); set(h, 'linewidth', 3);

axis([-1 1 -1 1 -1 1]); daspect([1 1 1]); box on;