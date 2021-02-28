function R = PQDprofile(P)
x = -30:30; y = -30:30;
[X, Y] = meshgrid(x, y);
R = angle(X+1i*Y);
R(find(R>(pi/3))) = 2*(pi/3)-R(find(R>(pi/3)));
R(find(R<0)) = -R(find(R<0));
R(find(R>(pi/3))) = -(2*(pi/3)-R(find(R>(pi/3))));
R(find(R<0)) = -R(find(R<0));
end
