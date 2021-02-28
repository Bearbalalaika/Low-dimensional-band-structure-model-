%ALL RESULTS APPEAR IN THE STRUCTURE Q
tic
YES = 1; NO = 0;
CB = YES;     CB_LEVELS = 2;
HH = NO;     HH_LEVELS = 2;
LH = NO;     LH_LEVELS = 2;

VB = YES;    VB_LEVELS = 26; VB_LEVELS = 10;

clear Q;

TETRAHEDER = YES;

SPHERICAL_BANDS = NO;
COUPLING_OFF = NO;
STRAIN = NO;
PIEZO = NO;

HKL = [1 1 1];

%%%DEFINE HIGH RESOLUTION COORDINATE SYSTEM ------------------------------
dx0 = 0.5; dy0 = 0.5; dz0 = 0.25;
x0 = -20:dx0:20;
y0 = -20:dy0:20;

r = 9;

if (r == 10.50);
    QWR_0p35 = 0.0525; QW_0p35 = 0.1500
    QWR_0p20 = 0.010;  QW_0p20 = 0.1010;
elseif (r == 8.00)
    QWR_0p35 = 0.0437; QW_0p35 = 0.1450
    QWR_0p20 = 0.0020; QW_0p20 = 0.0975;
end    


%OLD parameters
QWR_0p35 = 0.0458; QW_0p35 = 0.204;
QWR_0p20 = 0.0103; QW_0p20 = 0.075;

QWR_0p35 = 1.10*QWR_0p35;
QWR_0p20 = 1.00*QWR_0p20;

t1 = 7; %QD height
t2 = 7.5; %QD height

% conc = [0.30 0.30 0.30  32/dz0;   %Cladding
%         0.20 QW_0p20 QWR_0p20  t1/dz0;   %QD 1
%         0.30 0.30 0.30  32/dz0;]; %Cladding
% 
% 
% [C1, z0] = flat_pyramid_layers(x0, y0, conc, r);
 
 
% conc = [0.30 0.30 0.30  32/dz0;   %Cladding
%         0.20 QW_0p20 QWR_0p20  t2/dz0;   %QD 2
%         0.30 0.30 0.30  (32-(t2-t1))/dz0;]; %Cladding
 
 
% [C2, z0] = flat_pyramid_layers(x0, y0, conc, r);

% conc = [0.35 0.35 0.35  32/dz0;   %Cladding
%         0.20 QW_0p20 QWR_0p20  t1/dz0;   %QD 1
%        0.35 0.35 0.35  32/dz0;]; %Cladding
% 
% 
%[C1, z0] = flat_pyramid_layers(x0, y0, conc, r);
% 
%
% conc = [0.35 0.35 0.35  32/dz0;   %Cladding
%         0.20 QW_0p20 QWR_0p20  t2/dz0;   %QD 1
%         0.35 0.35 0.35  (32-(t2-t1))/dz0;]; %Cladding
% 
% 
% [C2, z0] = flat_pyramid_layers(x0, y0, conc, r);
%
 conc = [0.35 QW_0p35 QWR_0p35  32/dz0;   %Cladding
         0.20 QW_0p20 QWR_0p20  t1/dz0;   %QD 1
         0.35 QW_0p35 QWR_0p35  32/dz0;]; %Cladding
 
 
 [C1, z0] = flat_pyramid_layers(x0, y0, conc, r);
 
 
 conc = [0.35 QW_0p35 QWR_0p35  32/dz0;   %Cladding
         0.20 QW_0p20 QWR_0p20  t2/dz0;   %QD 1
         0.35 QW_0p35 QWR_0p35  (32-(t2-t1))/dz0;]; %Cladding
 
 
 [C2, z0] = flat_pyramid_layers(x0, y0, conc, r);


z0 = z0*dz0;
disp(max(z0));

C = C1;

    TP =tetra_profile(x0, y0);
    TP(find(TP<15)) = 15; TP = gaussBlur(TP, 5, 5, 5);
    TP = TP - min(TP(:));

    c0 = C(:,:,size(C,3));

    TPV = zeros(size(C));
    for i=1:size(C,3);
        TPV(:,:,i) = TP;
    end

    [X, Y, Z] = meshgrid(x0, y0, z0);
    C = interp3(X, Y, Z, C, X, Y, Z-TPV);

    for i=1:size(C,3);
        c_qd = C(:,:,i);    c_qd(isnan(c_qd)) = c0(isnan(c_qd));    C(:,:,i) = c_qd;
    end

    %C = gaussBlur(C, 2, 2, 2);
    clear TP X Y Z TPV c0 c_qd;

    %Diffusion model
    %Ld = 4/dz0; %diffusion length
    %C = gaussBlur(C, 0, 0, Ld);
C1 = C;

C = C2;

    TP =tetra_profile(x0, y0);
    TP(find(TP<15)) = 15; TP = gaussBlur(TP, 5, 5, 5);
    TP = TP - min(TP(:));

    c0 = C(:,:,size(C,3));

    TPV = zeros(size(C));
    for i=1:size(C,3);
        TPV(:,:,i) = TP;
    end

    [X, Y, Z] = meshgrid(x0, y0, z0);
    C = interp3(X, Y, Z, C, X, Y, Z-TPV);

    for i=1:size(C,3);
        c_qd = C(:,:,i);    c_qd(isnan(c_qd)) = c0(isnan(c_qd));    C(:,:,i) = c_qd;
    end

    %C = gaussBlur(C, 2, 2, 2);
    clear TP X Y Z TPV c0 c_qd;

    %Diffusion model
    %Ld = 4/dz0; %diffusion length
    %C = gaussBlur(C, 0, 0, Ld);
C2 = C;

mpfname = 'InAlGaAs_mparam_JAP89';

Q.hrpot.C_QD = single(0);
Q.hrpot.C_VQWR = single(0);
Q.hrpot.C_VQW = single(0);
Q.hrpot.C = single(0);
Q.hrpot.x = x0;
Q.hrpot.y = y0;
Q.hrpot.z = z0;
clear C



%%%DEFINE & LOAD MATERIAL PARAMETERS
if strcmp(mpfname, 'basic')
    disp('BASIC material parameters');
    mparam = generate_basic_mparam(0.065);
else
    mpfname = 'InAlGaAs_mparam_JAP89';
    disp(['Load material parameters: ' mpfname]);
    load(mpfname);
end
Q.mparam = mparam;
Q.mparam.file = mpfname; 
clear mpfname mparam;

%%%DEFINE Z-DIRECTION OF CRYSTAL, Y direction is [-h k 0] and X follow
Q.hkl = HKL;   %Use NaN for principal axes

if STRAIN
%%%Simulate strain relaxation, result: strain fields
Q.strain = call_strainND(Q, 1);

%%%Extract coordinate system for used for strain simulation
xs = Q.strain.x;
ys = Q.strain.y;
zs = Q.strain.z;
end

if PIEZO
%%%Calculate the charge density
Q.piezo.rho = piezo_charge(Q.strain, Q.mparam, xs, ys, zs, Q.hkl, 'INV');

%%%Calculate the piezoelectric potential (felt by electrons)
[Q.piezo.V, Q.piezo.Vext] = call_poissonND(xs, ys, zs, Q.strain.c1 + (1i)*Q.strain.c2, Q.piezo.rho, 1);  %Result in eV
end


%%%DEFINE LOW RESOLUTION COORDINATE SYSTEM ------------------------------

dx = 1; dy = 1; dz = 0.5;
x = -18:dx:18;
y = -18:dy:18;
z = -32:dz:32;


%Resample the potential on the new low-res grid
Q.lrpot.C = single(resampleND(x0, y0, z0, C2, x, y, z)); clear C2;
C1 = single(resampleND(x0, y0, z0, C1, x, y, z));
Q.lrpot.x = x;
Q.lrpot.y = y;
Q.lrpot.z = z;

Q.hrpot = 0;

Q.lrpot.C = Clone_structure_Z(Q.lrpot.C, -round(t1/dz+d), +round( (t1/dz+d)/2 ), C1); clear C1;


if VB
Q.VB.levels = VB_LEVELS;
Q.VB.spherical_approx = SPHERICAL_BANDS;
Q.VB.coupling_off = COUPLING_OFF;


if (isnan(Q.hkl)|(Q.hkl == [0 0 1]))
    Q.VB.hamiltonian = 'lutt3D_001'
elseif (Q.hkl == [1 1 1])
    Q.VB.hamiltonian = 'lutt3D_111';
else
    disp('WARNING: No Luttinger Hamiltonian defined for the choosen Z-axis!!');
end
[Q.VB.E, Q.VB.WF] = schrodND(Q, Q.VB.hamiltonian);
end

if CB
    Q.CB.levels = CB_LEVELS;
    Q.CB.hamiltonian = 'Hc'
    [Q.CB.E, Q.CB.WF] = schrodND(Q, Q.CB.hamiltonian);
end

if HH
    Q.HH.levels = HH_LEVELS;
    Q.HH.hamiltonian = 'HH3D_111'
    [Q.HH.E, Q.HH.WF] = schrodND(Q, Q.HH.hamiltonian);
end

if LH
    Q.LH.levels = LH_LEVELS;
    Q.LH.hamiltonian = 'LH3D_111'
    [Q.LH.E, Q.LH.WF] = schrodND(Q, Q.LH.hamiltonian);
end



toc
save(['UNEQUAL_' num2str(t1) 'nm_DQD_' num2str(d) 'dz_Barrier_NOVQWR'], 'Q')
