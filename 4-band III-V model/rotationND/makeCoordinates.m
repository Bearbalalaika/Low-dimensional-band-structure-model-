function s = makeCoordinates(s)

[Sy, Sx, Sz] = size(s.e_xx);

px = 1:(Sx+1); py = 1:(Sy+1); pz = 1:(Sz+1);
[P0x, P0y, P0z] = meshgrid(px, py, pz);  %Initialize xyz cooordinates

L= 1./(s.e_misfit+1);

Px = 0*P0x;
Py = 0*P0y;
Pz = 0*P0z;

x = 2; y = 1; z = 3;  %Dimensions

Px = Px - split_neighbours(s.e_xx, [-x y z])/2;
Py = Py - split_neighbours(s.e_yy, [ x -y z])/2;
Pz = Pz - split_neighbours(s.e_zz, [ x y -z])/2;

%Px = Px - split_neighbours(s.e_xy, [x -y z])/4;
%Py = Py - split_neighbours(s.e_xy, [-x y z])/4;

%Px = Px - split_neighbours(s.e_xz, [x y -z])/4;
%Pz = Pz - split_neighbours(s.e_xz, [-x y z])/4;

%Py = Py - split_neighbours(s.e_yz, [x y -z])/4;
%Pz = Pz - split_neighbours(s.e_yz, [x -y z])/4;

s.Px = Px+P0x;
s.Py = Py+P0y;
s.Pz = Pz+P0z;



P0 = mean_neighbours(diff(s.Px, 1, x), [y z]);
s.e_xx = 1-P0./L;

P0 = mean_neighbours(diff(s.Py, 1, y), [x z]);
s.e_yy = 1-P0./L;

P0 = mean_neighbours(diff(s.Pz, 1, z), [x y]);
s.e_zz = 1-P0./L;

P0 = mean_neighbours(diff(s.Px, 1, y), [x z]) + mean_neighbours(diff(s.Py, 1, x), [y z]);
s.e_xy = -0.5*P0./L;

P0 = mean_neighbours(diff(s.Px, 1, z), [x y]) + mean_neighbours(diff(s.Pz, 1, x), [y z]);
s.e_xz = -0.5*P0./L;

P0 = mean_neighbours(diff(s.Py, 1, z), [x y]) + mean_neighbours(diff(s.Pz, 1, y), [x z]);
s.e_yz = -0.5*P0./L;

