%out = schrodQW([0.25 0. 0.3],  20, [3 3 3], 0, 0);
%out = schrodQW(z, Al, levels, fields, kxy)
%
%width = 15; fields = 0; kxy = [0 0];
%z = -20:0.1:20;
%Al = zeros(size(z)); Al(find(abs(z)>(width/2))) = 0.33;
function out = schrodQW(z, Al, profile, levels, fields, kxy)
disp('Fields not implemented')
tic
if length(Al) == 1;
    dz = z(1);   
    width = Al;
    Al_QW = z(2);
    Al_barrier = z(3);
    
    z = -max(20, round((2.5*width))):dz:max(round((2.5*width)), 20);%spatial domain for calculation
    
    Al = Al_QW*ones(size(z)); 
    Al(find(z<(-width/2))) = Al_barrier; 
    Al(find(z>=(width/2))) = Al_barrier;
else
    dz = mean(diff(z));
end
%for rctangular QW profile=0 (or any other nr than 1,2,3)

%different potentil profiles

if profile==1 %triangular QW
    
    for i=0:1:width/dz
    Al(find(z==(-width/2+i*dz))) = Al_QW + i*dz/width*(Al_barrier-Al_QW);
    end
    Al(find(z<(-width/2))) = Al_barrier; 
    Al(find(z>=(width/2))) = Al_barrier;
    
end


if profile==2 %V-shape QW
    
    for i=0:1:width/(2*dz)
    Al(find(z==(-width/2+i*dz))) = Al_barrier - i*dz/(0.5*width)*(Al_barrier-Al_QW);
    end
    for i=0:1:width/(2*dz)
    Al(find(z==(i*dz))) = Al_QW + i*dz/(0.5*width)*(Al_barrier-Al_QW);
    end
    
    Al(find(z<(-width/2))) = Al_barrier; 
    Al(find(z>=(width/2))) = Al_barrier;
    
    Al=Al./(Al+8.6.*(1-Al)); %effective Al content (lateral confinement energy included)
end


if profile==3 %parabolic
    
    for i=0:1:width/dz+1
    Al(find(z==(-width/2+i*dz))) = Al_QW + (-width/2+i*dz)*(-width/2+i*dz)*4/(width*width)*(Al_barrier-Al_QW);
    end

    Al(find(z<(-width/2))) = Al_barrier; 
    Al(find(z>=(width/2))) = Al_barrier;
end


Eg = zeros(size(z));
    
for i=1:1:size(z)
    Eg(i) = 1.519 + 1.707.*Al(i) - 1.437.*Al(i).*Al(i) + 1.310.*Al(i).*Al(i).*Al(i);
end
    
out.Eg = Eg;

if nargin <3
   levels = [3 3 3];
elseif length(levels) == 1
   levels = [levels 0 0];
elseif length(levels) == 2
   levels = [levels 0];
end   

if nargin <5
    kx = 0;
    ky = 0;
else
    if length(kxy) == 1
        kx = kxy/sqrt(2);
        ky = kx;
    elseif length(kxy) == 2
        kx = kxy(1);
        ky = kxy(2);
    end
end

physconst
%load('c:\matlab\QWR\QWRs\AlGaAs_mparam2');
load('AlGaAs_mparam_JAP89'); 
%disp('Q = 0.7');  mparam.Q = 0.7;
%disp('Egap modified');  mparam.Egap = [0.22 1.36 1.5192];
%disp('Diamagnetic test material parameters'); load('e:\matlab\QWRs\AlGaAs_mparamDia');
%load('AlGaAs_mparam_PRB37_3130'); disp('PRB37_3130')

MASSE = polyval(mparam.me , Al);
MASSHH = polyval(mparam.mhh, Al);
MASSLH = polyval(mparam.mlh, Al);

POTF = 0;
POTMZ = hbar^2/(2*me*eV) * ((kx*1E9).^2 + (ky*1E9).^2);
POTME = POTMZ./MASSE;
POTMHH = POTMZ./MASSHH;
POTMLH = POTMZ./MASSLH;


POTE = mparam.Q*(polyval(mparam.Egap, Al)-mparam.Egap(length(mparam.Egap))) + POTF + POTME;
POTHH = (1-mparam.Q)*(polyval(mparam.Egap, Al)-mparam.Egap(length(mparam.Egap))) - POTF + POTMHH;
POTLH = (1-mparam.Q)*(polyval(mparam.Egap, Al)-mparam.Egap(length(mparam.Egap))) - POTF + POTMLH;

if levels(1)>0
	disp(['Electron levels ' num2str(levels(1)) '-------------------------']);
	[Ee, WFe] = schrod1D(z, POTE, MASSE, levels(1));
end

if levels(2)>0
	disp(['heavy hole levels ' num2str(levels(2)) '-------------------------']);
	[Ehh, WFhh] = schrod1D(z, POTHH, MASSHH, levels(2));
end
if levels(3)>0
	disp(['light hole levels ' num2str(levels(3)) '-------------------------']);
	[Elh, WFlh] = schrod1D(z, POTLH, MASSLH, levels(3));
end

comptime = toc;
%out.TEM = TEM;
out.mparam = mparam;

%out.cparam.B = B;
%out.cparam.F = F/100E3;
out.cparam.kxy = kxy;
%out.cparam.Al_QWR = Al_QWR;
%out.cparam.Al_VQW = Al_VQW;
%out.cparam.Al_barrier = Al_barrier;
out.cparam.Al = Al;
out.cparam.z = z;


out.cparam.Fpot = POTF;

    
out.cparam.scale = dz;
%out.cparam.reduction = factor;
%out.cparam.sizeOut = sizeOut;

%out.cparam.realreduction = sizeOut./size(POT);
%out.cparam.realscale = out.cparam.scale./out.cparam.realreduction;

%out.cparam.center = [data.parameters.center(1) data.parameters.center(2)];
out.cparam.comptime = comptime;

if levels(1) > 0
    out.E.E = Ee;
    out.E.WF = WFe;
end

if levels(2) > 0
    out.HH.E = Ehh;
    out.HH.WF = WFhh;
end

if levels(3) > 0
    out.LH.E = Elh;
    out.LH.WF = WFlh;
end
