function C = finite_VQWR(x, y, z, shape)
if nargin<3
    shape = x;
    x = -25:0.5:25;
    y = -25:0.5:25;
    z = -60:0.5:60;
end

%x = 0.30; VQW = x/(x+2.1*(1-x)); VQWR = x/(x+8.6*(1-x));
%1.518 = 1.36x + 0.22x^2

%OLD SERIES
%conc_bulk_top_bottom = 0.55;
%conc_bulk_center = 0.3;

%conc_VQW_top_bottom = concVQW(conc_bulk_top_bottom);
%conc_VQW_center = concVQW(conc_bulk_center);
%conc_VQWR_top_bottom = concVQWR(conc_bulk_top_bottom);
%conc_VQWR_center = concVQWR(conc_bulk_center);

%conc_VQWR_top_bottom = 0.09
%conc_VQW_center = 0.145
%conc_VQWR_center = 0.02

%NEW SERIES
conc_bulk_top_bottom = 0.35;
conc_bulk_center = 0.2;

conc_VQW_top_bottom = concVQW(conc_bulk_top_bottom);
conc_VQW_center = concVQW(conc_bulk_center);
conc_VQWR_top_bottom = concVQWR(conc_bulk_top_bottom);
conc_VQWR_center = concVQWR(conc_bulk_center);

conc_VQWR_top_bottom = 0.09*(0.35/0.55)*0.8
conc_VQW_center = 0.145*(0.2/0.3)*0.775
conc_VQWR_center = 0.02*(0.2/0.3)*0.775

t1 = shape;

r1 = 10;
r0 = 10;

d1 = 17;
d0 = 17;

r1 = 8;
r0 = r1;

d1 = 13;
d0 = d1;

%C = 0.30*buildSINGLEVQW(x, y, z, 12);
%C(find(C==0)) = 0.14;
%C = tcone(x, y, z, r1, 0, 0, -t1/2) - tcone(x, y, z, r1, 0, 0, +t1/2);


if 1
    T = tcone(x, y, z, r1, 0, 0, -t1/2) - tcone(x, y, z, r1, 0, 0, +t1/2);

    C = conc_bulk_top_bottom*ones(length(y), length(x), length(z));
    C(find(T==1)) =  conc_bulk_center;

    [C0, C1] = buildVQWs(x, y, z, d0, d1, 0, 0);
    C(find(C0&(~T))) = conc_VQW_top_bottom;
    C(find(C1&T)) = conc_VQW_center;

    [C0, C1] = buildVQWRs(x, y, z, r0, r1, t1, 0, 0, 0);
    C(find(C0&(~T))) = conc_VQWR_top_bottom;
    C(find(C1&T)) = conc_VQWR_center;
end

function x = concVQW(x)
x = x/(x+2.1*(1-x));

function x = concVQWR(x)
x = x/(x+8.6*(1-x));


function C = tcone(x, y, z, r1, x0, y0, z0)
[X, Y, Z] = meshgrid(x-x0, y-y0, z-z0);

%cone
%F = sqrt(2);
%C = zeros(size(X));
%R = F*sqrt(X.^2 + Y.^2);
%C(find((Z+F*r1)>R)) = 1;
%C(find(Z<0)) = 0;

%flat
C = ones(size(X));
C(find(Z<0)) = 0;

%strange
%FI = angle(X+1i*Y);%-pi*sign(Y+1E-5);
%p1 = 1+0.05; p2 = 0.20;
%R = R./(1-(1./(p1+1)).^p2+ (1./(p1+cos(3*FI))).^p2);
%C(find((Z+F*r1)>R)) = 1;
%C(find(Z<0)) = 0;



function [C0, C1] = buildVQWRs(x, y, z, r0, r1, t1, x0, y0, z0)
[X, Y, Z] = meshgrid(x-x0, y-y0, z-z0);
R = sqrt(X.^2 + Y.^2);
C0 = zeros(size(X));
C1 = zeros(size(X));
C1(find(R<r1)) = 1;
C0(find(R<r0)) = 1;



function [C0, C1] = buildVQWs(x, y, z, a0, a1, x0, y0)
C0 = buildVQW(x, y, z, a0, x0, y0);
C1 = buildVQW(x, y, z, a1, x0, y0);


function C = buildSINGLEVQW(x, y, z, a)
[X, Y, Z] = meshgrid(x, y, z);
L = a;
C = ones(size(X));
f = -0*pi/180;
M = sin(f)*X-cos(f)*Y;
C(find((abs(M)<(L/2)))) = 0;
