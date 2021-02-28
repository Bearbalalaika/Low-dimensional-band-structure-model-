%T.a = 10; T.type = 'tetra'; T.x0 = 0; T.y0 = 0; T.z0 = 0; C = tetra(T);
%figure; isosurface(C); daspect([1 1 1]); box on ;camlight left; lighting gouraud  
function [C, C2, C3] = tetra(x, y, z, shape)
if nargin<3
    shape = x;
    x = -20:0.5:20;
    y = -20:0.5:20;
    z = -20:0.5:20;
end
%
%
%z_sep = 10;
%C1B = buildtetra(x, y, z, a, 0, 0, 0-z_sep, 1, 0, 7);
%C2B = buildtetra(x, y, z, a, 0, 0, -2-z_sep, 1, 0, 10);


%C4 = buildtetra(x, y, z, a, -10, 10, 0, 1, 0);
%C5 = buildtetra(x, y, z, a, -10, -10, 0, 1, 0)

%C1B(find(C2B>0))=0;

%C = double(C1A|C1B);

if strcmp(shape.type, 'tetra')
    if ~isfield(shape, 'trunc')
        shape.trunc = inf;
    end
    C=buildtetra(x, y, z, shape.a, shape.x0, shape.y0, shape.z0, 1, 0, -shape.trunc);
    C2 = NaN; C3 = NaN;
end


%Fab APL
if strcmp(shape.type, 'fabprof')
    if ~isfield(shape, 'InQD')
        shape.InQD = 0;
    end
    h = shape.h;
    w = shape.w;
    t = h/8;
    %a_vqwr = 22;
    a_vqwr = 17;

    shift_x = -3;
    if h>0
        C1 = fabprof(w, h, t, x, y, -z, shift_x, 0, -12);
        if isfield(shape, 'noSQW')
            C1(:,:,find((z)>=0)) = 0;
        end
    else
        C1 = zeros(length(y), length(x), length(z));
    end    
    C2 = buildVQWR(x, y, z, sqrt(3/2)*a_vqwr, shift_x, 0);
    C2(find(C1>0)) = 0;
    C2 = C2*0.05;   % 5% Aluminium in VQW
    C3 = buildVQW(x, y, z, 0.5*sqrt(3/2)*a_vqwr, shift_x, 0);
    C3(find(C1>0)) = 0;
    C3 = 0.18*C3;  % 18% Aluminium in VQW
    C = C1;
    C(find(C1>0)) =  shape.InQD*(1i);   % ?% Indium in QD
    C(find(C1==0)) = 0.30;   % 30% Aluminium in barrier
    %C = (1-C)*0.3;  % 0% In in QD, 30% Al in barrier        
    
    %if isfield(shape, 'blur') C = gaussBlur(C, shape.blur, shape.blur, shape.blur); end;
    %figure; plot(reshape(C(41,41,:), 81, 1))
end


if strcmp(shape.type, 'dalessi')
    a = shape.a;
    C2 = buildVQWR(x, y, z, a, -5, 0);
    
    C = 0.3 + 0*C2;

    C2 = C2*0.06;
    C3 = buildVQW(x, y, z, a/2, -5, 0);
    C3 = 0.2*C3;
    %if isfield(shape, 'blur') C = gaussBlur(C, shape.blur, shape.blur, shape.blur); end;
    %figure; plot(reshape(C(41,41,:), 81, 1))
end



function C = cut(x, y, z, a, x0, y0, fi)
if nargin<5 x0 = 0; end
if nargin<6 y0 = 0; end
if nargin<7 fi = 60; end
fi = fi*pi/180;
[X, Y] = meshgrid(x-x0, y-y0);
C0 = ones(size(X));
C = zeros(length(y), length(x), length(z));

Cx = -1/cos(fi);
Cy = 1/sin(fi);
P1 = Cx*X+Cy*Y;
C0(find((Y>0)&(P1>0))) = 0;
for i=1:length(z)
	C(:,:,i) = C0;
end


function C=buildVQWR(x, y, z, a, x0, y0, Vin, Vout)
if nargin<5 x0 = 0; end
if nargin<6 y0 = 0; end
if nargin<8 Vin = 1; end
if nargin<9 Vout = 0; end
C0 = buildtetra(x, y, 0, a, x0, y0,0, Vin, Vout, 1, -1);
C = zeros(length(y), length(x), length(z));
for i=1:length(z)
	C(:,:,i) = C0;
end

function C = buildVQW(x, y, z, a, x0, y0)
[X, Y, Z] = meshgrid(x-x0, y-y0, z);
L = a;
C = ones(size(X));
f = -0*pi/180;
M = sin(f)*X-cos(f)*Y;
C(find((X>=0)&(abs(M)<(L/2)))) = 0;

f = -60*pi/180;
M = sin(f)*X-cos(f)*Y;
C(find((Y>=0)&(abs(M)<(L/2)))) = 0;

f = 60*pi/180;
M = sin(f)*X-cos(f)*Y;
C(find((Y<=0)&(abs(M)<(L/2)))) = 0;
C = 1-C;


function C=fabprof(w, h, t, x, y, z, x0, y0, z0)
if nargin<7 x0 = 0; end
if nargin<8 y0 = 0; end
if nargin<9 z0 = 0; end
sin54_73 = 0.8165;
w0 = w;
a = 40;
trunc_top1 = 2*sqrt(2)/3*w0; w1 = w0+3/2*(h/sqrt(2)-t/sin54_73);
trunc_top2 = 2*sqrt(2)/3*w1;
dtrunc = trunc_top2-trunc_top1;
C = buildtetra(x, y, z, a, x0, y0, z0, 1, 0, trunc_top1);
C2 = buildtetra(x, y, z, a, x0, y0, z0-h+dtrunc, 1, 0, trunc_top1+dtrunc);
C(find(C2>0))=0;

function C=buildtetra(x, y, z, a, x0, y0, z0, Vin, Vout, ttop, tbottom)
if nargin<5 x0 = 0; end
if nargin<6 y0 = 0; end
if nargin<7 z0 = 0; end
if nargin<8 Vin = 1; end
if nargin<9 Vout = 0; end
if nargin<10 ttop = -inf; end
if nargin<11 tbottom = 0; end

[X, Y, Z] = meshgrid(x-x0, y-y0, z-z0);

d = sqrt(6)/3*a;

C1x = d/(sqrt(3)*a/3);
C1y = d/(a/3);
C1z = d/(sqrt(6)*a/3);
P1 = C1x*X + C1y*Y + C1z*Z -d;

C2x = d/(sqrt(3)*a/3);
C2y = -d/(a/3);
C2z = d/(sqrt(6)*a/3);
P2 = C2x*X + C2y*Y + C2z*Z -d;


C3x = d/(-sqrt(3)*a/6);
C3y = 0;
C3z = d/(sqrt(6)*a/3);
P3 = C3x*X + C3y*Y + C3z*Z - d;

C = ones(size(X))*Vout;
%C(find((P1<0)&(P2<0)&(P3<0)&(Z>=tbottom)&(Z<(d-ttop)))) = Vin;
C(find((P1<=0)&(P2<=0)&(P3<=0)&(Z>=tbottom)&(Z<=(d-ttop)))) = Vin; disp('equal')


