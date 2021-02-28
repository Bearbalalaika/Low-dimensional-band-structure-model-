%function H = schrod3D
%Comparasion with finitediff gives identical results for x-y if
%adoint of conc is used: for separate operatos Dx2 Dy2 Dxy Dxapprx Dyapprx
%All DIFF combination verified with the TESTMODULE for dxyz = 1;

function [E, WF] = schrodND(inp, hamiltonian);
conc = inp.lrpot.C;
x = inp.lrpot.x;
y = inp.lrpot.y;
z = inp.lrpot.z;

kx = 0;     %MUST WRITE HANLDES FOR THIS
ky = 0;     %SO THEY CAN BE INCLUDED IN THE INPUT STRUCT
kz = 0;
y_sym = 0

BC = [];

physconst;

Stot = length(x)*length(y)*length(z);
Sx = length(x); Sy = length(y); Sz = length(z);

[X, Y, Z] = meshgrid(x, y, z);

dim = ~~(Sx-1) + ~~(Sy-1) + ~~(Sz-1);
%dim = 3;

disp(['Dimensions: ' num2str(dim)]);

CB = 0; HH = 0; LH = 0;
if strcmp(hamiltonian, 'Hc') [H, n] = H_CB(inp, X, Y, Z); CB = 1; end
if strcmp(hamiltonian, 'HH') [H, n] = H_HH(inp, X, Y, Z); HH = 1; end
if strcmp(hamiltonian, 'LH') [H, n] = H_LH(inp, X, Y, Z); LH = 1; end
if strcmp(hamiltonian, 'lutt1D_001')   [H, n] = lutt1D_001(inp, X, Y, Z, k1, k2); end
if strcmp(hamiltonian, 'lutt2Dky_001') [H, n] = lutt2D_ky_001(inp, X, Y, Z); dim = 3; end
if strcmp(hamiltonian, 'lutt2Dkz_001') [H, n] = lutt2D_kz_001(inp, X, Y, Z); end
if strcmp(hamiltonian, 'lutt3D_001')   [H, n] = lutt3D_001(inp, X, Y, Z); end

if strcmp(hamiltonian, 'lutt1D_111')   [H, n] = lutt1D_111(inp, X, Y, Z, k1, k2); end
if strcmp(hamiltonian, 'lutt3D_111')   [H, n] = lutt3D_111(inp, X, Y, Z); end
if strcmp(hamiltonian, 'lutt2Dkz_111') [H, n] = lutt2D_kz_111(inp, X, Y, Z); end 

if strcmp(hamiltonian, 'HH3D_111') [H, n] = HH3D_111(inp, X, Y, Z); HH = 1; end
if strcmp(hamiltonian, 'LH3D_111') [H, n] = LH3D_111(inp, X, Y, Z); LH = 1; end

if CB
    levels = inp.CB.levels;
    if isfield(inp.lrpot, 'HCfname')
        HC = load(inp.lrpot.HCfname);
        H = H + HC.M; clear HC
    end
    if isfield(inp.CB, 'MAG')
        HM = load(inp.CB.MAG);
        H = H + HM.M; clear HM
    end
elseif HH
    disp('No Bir-Pikus or Magnetic Hamiltonian loaded for HH')
    levels = inp.HH.levels;
elseif LH
    disp('No Bir-Pikus or Magnetic Hamiltonian loaded for LH')
    levels = inp.LH.levels;    
else
    levels = inp.VB.levels;
    if isfield(inp.lrpot, 'BPfname')
        BP = load(inp.lrpot.BPfname);
        H = H + BP.M; clear BP
    end
    if isfield(inp.VB, 'MAG')
        HM = load(inp.VB.MAG);
        H = H + HM.M; clear HM
    end
end

%Clear everything not nesessary
clear X Y Z;
clear inp;
disp('Find eigenvalues');


OPTIONS.disp = 1;
OPTIONS.tol = 1E-25; %-11
OPTIONS.maxit = 2500;
disp('Save sparse Hamiltonian -> schrodND_H'); save schrodND_H H;
[WF, E] = eigs(H/meV, levels, 'SR', OPTIONS); E = diag(E);
clear H;

[E, I] = sort(E);
 
[X, Y, Z] = meshgrid(x, y, z);
iLx = invLx(X, Y, Z); iLx = kron(iLx,diag(ones(1,n), 0));
iLy = invLy(X, Y, Z); iLy = kron(iLy,diag(ones(1,n), 0));
iLz = invLz(X, Y, Z); iLz = kron(iLz,diag(ones(1,n), 0));
WF = single(iLz*iLy*iLx*WF);

clear iLx iLy iLz X Y Z;
disp('Save unreshaped wavefunction -> schrodND_result'); save schrodND_result WF E;

disp('Reshape wavefunctions')
if (dim==3)
    outWF = zeros(Sy, Sx, Sz*n, length(I), 'single');
    WF0   = zeros(Sy, Sx, Sz*n, 'single');
end

if (dim==2)
    outWF = zeros(Sy, Sx*n, length(I), 'single');
    WF0   = zeros(Sy, Sx*n, 'single');
end

if (dim==1)
    outWF = zeros(Sy*n, length(I), 'single');
    WF0   = zeros(Sy*n, 'single');
end

for i=1:length(I)
    
    if (dim==1) WF0 = reshape(WF(:,I(i)), Sy*n, 1);          outWF(:,i) = WF0; end
    if (dim==2) WF0 = reshape(WF(:,I(i)), Sy, Sx*n);     outWF(:,:,i) = WF0; end
    if (dim==3) WF0 = reshape(WF(:,I(i)), Sy, Sx, Sz*n); outWF(:,:,:,i) = WF0; end
end

clear WF0;
WF = outWF;

disp('Interpolation for irregular mesh: Switched off!')
if 0
dx = 5; dy = 5; dz =5;
x0 = min(x):dx:max(x); 
y0 = min(y):dy:max(y); 
z0 = min(z):dz:max(z); 
WF = interpWF(WF, X, Y, Z, x0, y0, z0, n, E);
end            

function WF = interpWF(WF, X, Y, Z, x0, y0, z0, n, E)
[Sy, Sx, Sz] = size(X); Sy0 = length(y0); Sz0 = length(z0);
dim = ~~(Sx-1) + ~~(Sy-1) + ~~(Sz-1);
Sx0 = length(x0); 
disp('Interpolate wavefunctions')
outWF0 = [];
WF_new = [];
for i=1:length(E)
    outWF = zeros(Sy0, n*Sx0, Sz0);
    for j=1:n
        if (dim==1)
            A0=WF(:,i); A = A0((1+Sx*(j-1)):(j*Sx));
            A = interp1(x, A, x0)';
            outWF((1+Sx0*(j-1)):(j*Sx0)) =  A;
        end
        if (dim==2)
            A0=WF(:,:,i); A = A0(:, (1+Sx*(j-1)):(j*Sx));
            [X0, Y0] = meshgrid(x0, y0);
            A = interp2(X(:,:,1), Y(:,:,1), A, X0, Y0);
            outWF(:,(1+Sx0*(j-1)):(j*Sx0)) = A;
        end
        if (dim==3)
            A0=WF(:,:,:,i); A = A0(:, (1+Sx*(j-1)):(j*Sx),:);
            [X0, Y0, Z0] = meshgrid(x0, y0, z0);
            A = interp3(X(:,:,:), Y(:,:,:),Z(:,:,:), A, X0, Y0, Z0);
            outWF(:,(1+Sx0*(j-1)):(j*Sx0),:) = A;
        end
    end
    outWF = outWF/sqrt(sum(sum(sum(outWF.^2))));
    if (dim==1) WF_new(:,i) = outWF; end
    if (dim==2) WF_new(:,:,i) = outWF; end
    if (dim==3) WF_new(:,:,:,i) = outWF; end

end
WF = WF_new;

function [M, n] = H_CB(inp, X, Y, Z)
    disp('Build cond. band Hamiltonian')
    conc = inp.lrpot.C;
    mparam = inp.mparam;
    
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    
    mass = real(polyval(mparam.me, conc));
    
    V = bandAlignment(conc, mparam, 'CB'); 
            
    piezo = 0; Epot = 0;
    if(isfield(inp, 'lrpiezo'))
        disp('Piezofield found: Add piezo potential');
        piezo = inp.lrpiezo.V;
    end
    
    if(isfield(inp.CB, 'Epot'))
        if inp.CB.Epot disp('Electric field found: Add electric potential'); end
        Epot = inp.CB.Epot;
    end
    
    [M, n] = SB(X, Y, Z, V, mass, Epot + piezo);
       
function [M, n] = H_HH(inp, X, Y, Z)
    disp('Build Heavy-hole Hamiltonian')
    conc = inp.lrpot.C;
    mparam = inp.mparam;
    
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    
    mass = real(polyval(mparam.mhh, conc));
    
    V = bandAlignment(conc, mparam, 'VB');

    piezo = 0; Epot = 0;
    if(isfield(inp, 'lrpiezo'))
        disp('Piezofield found: Add piezo potential');
        piezo = -inp.lrpiezo.V;
    end
    
    if(isfield(inp.HH, 'Epot'))
        if inp.HH.Epot disp('Electric field found: Add electric potential'); end
        Epot = inp.HH.Epot;
    end
    
    
    [M, n] = SB(X, Y, Z, V, mass, Epot + piezo);
    

    function [M, n] = H_LH(inp, X, Y, Z)
    disp('Build light-hole Hamiltonian')
    conc = inp.lrpot.C;
    mparam = inp.mparam;
    
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    
    mass = real(polyval(mparam.mlh, conc));
    
    V = bandAlignment(conc, mparam, 'VB');
    
    piezo = 0; Epot = 0;
    if(isfield(inp, 'lrpiezo'))
        disp('Piezofield found: Add piezo potential');
        piezo = -inp.lrpiezo.V;
    end
    
    if(isfield(inp.LH, 'Epot'))
        if inp.LH.Epot disp('Electric field found: Add electric potential'); end
        Epot = inp.LH.Epot;
    end
    
    [M, n] = SB(X, Y, Z, V, mass, Epot + piezo);
    
function [M, n] = SB(X, Y, Z, V, mass, Epot)              %Singleband hamiltonian
    disp('Single band')
    [Sy, Sx, Sz] = size(V);
    Stot = Sx*Sy*Sz;
    
    physconst;
 
    V = V*eV;
    if(nargin>5)
        disp('Electric potential found: Add electric potential potential');
        V = V + Epot*eV;
    end
 
    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
 
    
    C = -hbar^2/(2*me*(1E-9)^2);
    M = spalloc(Stot, Stot, 7*Stot);
    if Sx>1; M = M + D2x(1./mass, X, Y, Z); end
    if Sy>1; M = M + D2y(1./mass, X, Y, Z); end
    if Sz>1; M = M + D2z(1./mass, X, Y, Z); end
    M = C*M + V;
    n = 1;
    
    
function [M, n] = lutt1D_001(conc, X, Y, Z, k1, k2)
    disp('Build 1D Lutt. Hamiltonian (growth:Z[001] X[100] Y[010]) (approximate single & mixed derivatives)')
    disp('See PRB 53 1507 (1996))')
      
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    mparam = inp.mparam; physconst
    
    [lp1, lp2, lp3] = luttingerParameters(conc, mparam, 'INV');
    LP1 = spdiags(reshape(lp1, Stot, 1), 0, Stot, Stot);
    LP2 = spdiags(reshape(lp2, Stot, 1), 0, Stot, Stot);
    LP3 = spdiags(reshape(lp3, Stot, 1), 0, Stot, Stot);
      
    V = bandAlignment(conc, mparam, 'VB')*eV;

    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
    
    if (length(k1)>1) k1 = spdiags(reshape(k1   , Stot, 1), 0, Stot, Stot); end
    if (length(k2)>1) k2 = spdiags(reshape(k2   , Stot, 1), 0, Stot, Stot); end
    
    C = -hbar^2/(2*me*(1E-9)^2);
 
    HH =    - k1.^2*(LP1+LP2);
    HH = HH - k2.^2.*(LP1+LP2);
    HH = HH + D2z(lp1-2*lp2, X, Y, Z) +V/C;
    
    LH =    - k1.^2*(LP1-LP2);
    LH = LH - k2.^2*(LP1-LP2);
    LH = LH + D2z(lp1+2*lp2, X, Y, Z) +V/C;
    
    f1 = (lp3).*1i*(k1-1i*k2); F1 = spdiags(reshape(f1   , Stot, 1), 0, Stot, Stot);
    f2 = (lp2+lp3)/2.*(k1-1i*k2).^2; F2 = spdiags(reshape(f2   , Stot, 1), 0, Stot, Stot);
    f3 = (lp2-lp3)/2.*(k1+1i*k2).^2; F3 = spdiags(reshape(f3   , Stot, 1), 0, Stot, Stot);
    
    %F3 = F3*0; %axial approx;
    
    R =     2i*sqrt(3)*(F1*Dzapprx(X, Y, Z) + Dzapprx(X, Y, Z)*F1)/2;
    
    S =     -sqrt(3)*(F2 + F3);
    
    %B = -S  C = R
    O = sparse(length(HH), length(HH));
    M = C*[HH   R  S   O;
            R' LH  O   S;
            S'  O  LH -R;
            O   S' -R' HH];   
    n = 4;

    
function [M, n] = lutt1D_111(conc, X, Y, Z, k1, k2)
    disp('Build 1D Lutt. Hamiltonian (growth:Z[111] X[11-2] Y[1-10]) (approximate single & mixed derivatives)')
    disp('See Semicond. Sci. Technol. 4 904-909 (1989)')
    
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    mparam = inp.mparam; physconst;
    
    [lp1, lp2, lp3] = luttingerParameters(conc, mparam, 'INV');   
    
    LP1 = spdiags(reshape(lp1, Stot, 1), 0, Stot, Stot);
    LP2 = spdiags(reshape(lp2, Stot, 1), 0, Stot, Stot);
    LP3 = spdiags(reshape(lp3, Stot, 1), 0, Stot, Stot);
    
    V = bandAlignment(conc, mparam, 'VB')*eV;
    
    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
    
    if (length(k1)>1) k1 = spdiags(reshape(k1   , Stot, 1), 0, Stot, Stot); end
    if (length(k2)>1) k2 = spdiags(reshape(k2   , Stot, 1), 0, Stot, Stot); end
    
    C = -hbar^2/(2*me*(1E-9)^2);
 
    HH =    - k1.^2*(LP1+LP3);
    HH = HH - k2.^2.*(LP1+LP3);
    HH = HH + D2z(lp1-2*lp3, X, Y, Z) +V/C;
    
    LH =    - k1.^2*(LP1-LP3);
    LH = LH - k2.^2*(LP1-LP3);
    LH = LH + D2z(lp1+2*lp3, X, Y, Z) +V/C;
    
    f1 = (2*lp2+lp3).*(k1-1i*k2); F1 = spdiags(reshape(f1   , Stot, 1), 0, Stot, Stot);
    f2 = (lp2-lp3).*(k1+1i*k2).^2; F2 = spdiags(reshape(f2   , Stot, 1), 0, Stot, Stot);
    f3 = (lp2+2*lp3).*(k1-1i*k2).^2; F3 = spdiags(reshape(f3   , Stot, 1), 0, Stot, Stot);
    f4 = (lp2-lp3).*(k1+1i*k2); F4 = spdiags(reshape(f4   , Stot, 1), 0, Stot, Stot);
    
    
    B =     -2/sqrt(3)*(F1*Dzapprx(X, Y, Z) + Dzapprx(X, Y, Z)*F1)/2;
    B = B + 2*i/sqrt(6)*F2;
    
    R =   + 1/sqrt(3)*F3;
    R = R + 4i/sqrt(6)*(F4*Dzapprx(X, Y, Z) + Dzapprx(X, Y, Z)*F4)/2;
    
    %B = -S  C = R
    O = sparse(length(HH), length(HH));
    M = C*[HH   B  R   O;
            B' LH  O   R;
            R'  O  LH -B;
            O   R' -B' HH];   
    n = 4;

function [M, n] = lutt3D_001(inp, X, Y, Z)
    disp('Build 3D Lutt. Hamiltonian (growth:Z[001] X[100] Y[010]) (approximate single & mixed derivatives)')
    disp('G. Fishman PRB 52 11132 (1995) / U. Ekenberg Phys. Rev B. 32 5138 (1985) / T. Ikonic Phys. Rev. B 46 (1992)')
    %Note Fishman uses X[1 1 0] Y[1 -1 0]! Use Ekenberg or Ikonic!!! 060915

    conc = inp.lrpot.C;

    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    mparam = inp.mparam; physconst
    
    [lp1, lp2, lp3] = luttingerParameters(conc, mparam, 'INV');
    
    if isfield(inp.VB, 'spherical_approx')
        if inp.VB.spherical_approx  lp3 = (lp2+lp3)/2; lp2 = lp3; disp('----> Spherical band approximaion'); end
    end
    
    V = bandAlignment(conc, mparam, 'VB')*eV;
    
    if(isfield(inp.VB, 'Epot'))
       if inp.VB.Epot disp('Electric field found: Add electric potential'); end
       V = V + inp.VB.Epot*eV;
    end
    
    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
    C = -hbar^2/(2*me*(1E-9)^2);

    HH =      D2x(lp1+lp2, X, Y, Z);
    HH = HH + D2y(lp1+lp2, X, Y, Z);
    HH = HH + D2z(lp1-2*lp2, X, Y, Z) + V/C;
    
    LH =      D2x(lp1-lp2, X, Y, Z);
    LH = LH + D2y(lp1-lp2, X, Y, Z);
    LH = LH + D2z(lp1+2*lp2, X, Y, Z) + V/C;
    
    S =   + 2*sqrt(3)*(       (Dxz(lp3, X, Y, Z)+Dzx(lp3, X, Y, Z))/2);
    S = S + 2*sqrt(3)*( (-1i)*(Dyz(lp3, X, Y, Z)+Dzy(lp3, X, Y, Z))/2);
    
    tf1 = 1; tf2 = 1;
    tf1 = -1; disp('Transform G. Fishman -> U. Ekenberg');
    %tf2 = -1; disp('Transform G. Fishman -> Z. Ikonic');

    R =   - sqrt(3)*(          D2x( (lp3+lp2)/2, X, Y, Z) - D2y( (lp3+lp2)/2, X, Y, Z)) *tf2;
    R = R - sqrt(3)*(   (-1i)*(Dxy( (lp3+lp2)/2, X, Y, Z) + Dyx( (lp3+lp2)/2, X, Y, Z)))*tf2;
    R = R - sqrt(3)*(          D2x( (lp3-lp2)/2, X, Y, Z) - D2y( (lp3-lp2)/2, X, Y, Z)) *tf1;
    R = R - sqrt(3)*(   (+1i)*(Dxy( (lp3-lp2)/2, X, Y, Z) + Dyx( (lp3-lp2)/2, X, Y, Z)))*tf1;
    
    O = sparse(length(HH), length(HH));
    
    if isfield(inp.VB, 'coupling_off')
        if inp.VB.coupling_off  S = O; R = O; disp('----> Valence band coupling switched off'); end
    end
    
    M = C*[HH     -S     R      O;
           -S'    LH     O      R;
            R'     O    LH      S;
            O      R'    S'    HH];   
    n = 4;


    function [M, n] = lutt2D_ky_001(inp, X, Y, Z)
    disp('Build 2D Lutt. Hamiltonian (growth:Z[001] X[100] Y[010]:ky) (approximate single & mixed derivatives)')
    disp('G. Fishman PRB 52 11132 (1995) / U. Ekenberg Phys. Rev B. 32 5138 (1985) / T. Ikonic Phys. Rev. B 46 (1992)')
    %Note Fishman uses X[1 1 0] Y[1 -1 0]! Use Ekenberg or Ikonic!!! 060915
    conc = inp.lrpot.C;

    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    mparam = inp.mparam; physconst
    
    [lp1, lp2, lp3] = luttingerParameters(conc, mparam, 'LIN');
    LP1 = spdiags(reshape(lp1, Stot, 1), 0, Stot, Stot);
    LP2 = spdiags(reshape(lp2, Stot, 1), 0, Stot, Stot);
    LP3 = spdiags(reshape(lp3, Stot, 1), 0, Stot, Stot);
     
    if isfield(inp.VB, 'spherical_approx')
        if inp.VB.spherical_approx  lp3 = (lp2+lp3)/2; lp2 = lp3; disp('----> Spherical band approximaion'); end
    end

    %d/dz = i*ky
    if isfield(inp.VB, 'ky') iky = (1i)*inp.VB.ky; disp(['kz = ' num2str(iky/(1i)) ' nm-1']); else ikz = 0; end
       
    
    V = bandAlignment(conc, mparam, 'VB')*eV;
    
    %[max(max(max(V))) min(min(min(V)))]/eV
    %slut
    
    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
    C = -hbar^2/(2*me*(1E-9)^2);
    

    HH =      D2x(lp1+lp2, X, Y, Z);
    HH = HH + iky*iky*(LP1+LP2);
    HH = HH + D2z(lp1-2*lp2, X, Y, Z) + V/C;
    
    LH =      D2x(lp1-lp2, X, Y, Z);
    LH = LH + iky*iky*(LP1-LP2);
    LH = LH + D2z(lp1+2*lp2, X, Y, Z) + V/C;
    
    tf1 = 1; tf2 = 1;
    %tf1 = -1; disp('Transform G. Fishman -> U. Ekenberg');
    %tf2 = -1; disp('Transform G. Fishman -> Z. Ikonic');

    
    S =   + 2*sqrt(3)*(       (Dxz(lp3, X, Y, Z)   + Dzx(lp3, X, Y, Z))/2);
    S = S + 2*sqrt(3)*( (-1i)*(Dz(X, Y, Z)*LP3*iky + iky*LP3*Dz(X, Y, Z))/2);
    
    
    R =   - sqrt(3)*(          D2x( (lp3+lp2)/2, X, Y, Z)  - iky*iky*(LP3+LP2)/2         )*tf2;
    R = R - sqrt(3)*(   (-1i)*(Dx(X, Y, Z)*(LP3+LP2)*iky/2 + iky*(LP3+LP2)*Dx(X, Y, Z)/2))*tf2;
    R = R - sqrt(3)*(          D2x( (lp3-lp2)/2, X, Y, Z)  - iky*iky*(LP3-LP2)/2         )*tf1;
    R = R - sqrt(3)*(   (+1i)*(Dx(X, Y, Z)*(LP3-LP2)*iky/2 + iky*(LP3-LP2)*Dx(X, Y, Z)/2))*tf1; 
    
    O = sparse(length(HH), length(HH));
    
    if isfield(inp.VB, 'coupling_off')
        if inp.VB.coupling_off  S = O; R = O; disp('----> Valence band coupling switched off'); end
    end
    

        
    M = C*[HH     -S    R      O;
           -S'   LH     O      R;
            R'      O    LH      S;
            O       R'   S'     HH];
        
    n = 4;

    
function [M, n] = lutt2D_kz_001(inp, X, Y, Z)
    disp('Build 2D Lutt. Hamiltonian (growth:Z[001]:kz X[100] Y[010]) (approximate single & mixed derivatives)')
    disp('G. Fishman PRB 52 11132 (1995) / U. Ekenberg Phys. Rev B. 32 5138 (1985) / T. Ikonic Phys. Rev. B 46 (1992)')
    %Note Fishman uses X[1 1 0] Y[1 -1 0]! Use Ekenberg or Ikonic!!! 060915
    
    conc = inp.lrpot.C;

    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    mparam = inp.mparam; physconst
    
    [lp1, lp2, lp3] = luttingerParameters(conc, mparam, 'INV');
    LP1 = spdiags(reshape(lp1, Stot, 1), 0, Stot, Stot);
    LP2 = spdiags(reshape(lp2, Stot, 1), 0, Stot, Stot);
    LP3 = spdiags(reshape(lp3, Stot, 1), 0, Stot, Stot);
     
    if isfield(inp.VB, 'spherical_approx')
        if inp.VB.spherical_approx  lp3 = (lp2+lp3)/2; lp2 = lp3; disp('----> Spherical band approximaion'); end
    end

    %d/dz = i*kz
    if isfield(inp.VB, 'kz') ikz = (1i)*inp.VB.kz; disp(['kz = ' num2str(ikz/(1i)) ' nm-1']); else ikz = 0; end
       
    
    V = bandAlignment(conc, mparam, 'VB')*eV;
    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
    C = -hbar^2/(2*me*(1E-9)^2);
    
    HH =      D2x(lp1+lp2, X, Y, Z);
    HH = HH + D2y(lp1+lp2, X, Y, Z); 
    HH = HH + ikz*ikz*(LP1-2*LP2) + V/C;
    
    LH =      D2x(lp1-lp2, X, Y, Z);
    LH = LH + D2y(lp1-lp2, X, Y, Z);
    LH = LH + ikz*ikz*(LP1+2*LP2) + V/C;
        
    tf1 = 1; tf2 = 1;
    %tf1 = -1; disp('Transform G. Fishman -> U. Ekenberg');
    %tf2 = -1; disp('Transform G. Fishman -> Z. Ikonic');

    S =   + 2*sqrt(3)*(       (Dx(X, Y, Z)*LP3*ikz + ikz*LP3*Dx(X, Y, Z))/2);
    S = S + 2*sqrt(3)*( (-1i)*(Dy(X, Y, Z)*LP3*ikz + ikz*LP3*Dy(X, Y, Z))/2);
    
    R =   - sqrt(3)*(          D2x( (lp3+lp2)/2, X, Y, Z) - D2y( (lp3+lp2)/2, X, Y, Z)) *tf2;
    R = R - sqrt(3)*(   (-1i)*(Dxy( (lp3+lp2)/2, X, Y, Z) + Dyx( (lp3+lp2)/2, X, Y, Z)))*tf2;
    R = R - sqrt(3)*(          D2x( (lp3-lp2)/2, X, Y, Z) - D2y( (lp3-lp2)/2, X, Y, Z)) *tf1;
    R = R - sqrt(3)*(   (+1i)*(Dxy( (lp3-lp2)/2, X, Y, Z) + Dyx( (lp3-lp2)/2, X, Y, Z)))*tf1;

    O = sparse(length(HH), length(HH));
    
    if isfield(inp.VB, 'coupling_off')
        if inp.VB.coupling_off  S = O; R = O; disp('----> Valence band coupling switched off'); end
    end
    

        
    M = C*[HH     -S    R      O;
           -S'   LH     O      R;
            R'      O    LH      S;
            O       R'   S'     HH];
        
    n = 4;

function [M, n] = HH3D_111(inp, X, Y, Z)
    disp('Build 3D HH Hamiltonian (growth:Z[1 1 1] X[1 1 -2] Y[1 -1 0]) (approximate single & mixed derivatives)')

    conc = inp.lrpot.C;

    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    mparam = inp.mparam; physconst
    
    [lp1, lp2, lp3] = luttingerParameters(conc, mparam, 'INV');
       
    %if isfield(inp.VB, 'spherical_approx')
    %    if inp.VB.spherical_approx  lp3 = (lp2+lp3)/2; lp2 = lp3; disp('----> Spherical band approximaion'); end
    %end
    
    V = bandAlignment(conc, mparam, 'VB')*eV;
    if(isfield(inp.HH, 'Epot'))
        if inp.HH.Epot disp('Electric field found: Add electric potential'); end
        V = V + inp.HH.Epot*eV;
    end

    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
    C = -hbar^2/(2*me*(1E-9)^2);

    HH =      D2x(lp1+lp3, X, Y, Z);
    HH = HH + D2y(lp1+lp3, X, Y, Z);
    HH = HH + D2z(lp1-2*lp3, X, Y, Z) +V/C;
    
    M = C*HH;   
    n = 1;    

function [M, n] = LH3D_111(inp, X, Y, Z)
    disp('Build 3D HH Hamiltonian (growth:Z[1 1 1] X[1 1 -2] Y[1 -1 0]) (approximate single & mixed derivatives)')

    conc = inp.lrpot.C;

    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    mparam = inp.mparam; physconst
    
    [lp1, lp2, lp3] = luttingerParameters(conc, mparam, 'INV');
       
    %if isfield(inp.VB, 'spherical_approx')
    %    if inp.VB.spherical_approx  lp3 = (lp2+lp3)/2; lp2 = lp3; disp('----> Spherical band approximaion'); end
    %end
    
    V = bandAlignment(conc, mparam, 'VB')*eV;
    
    if(isfield(inp.LH, 'Epot'))
        if inp.LH.Epot disp('Electric field found: Add electric potential'); end
        V = V + inp.LH.Epot*eV;
    end

    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
    C = -hbar^2/(2*me*(1E-9)^2);

    LH =      D2x(lp1-lp3, X, Y, Z);
    LH = LH + D2y(lp1-lp3, X, Y, Z);
    LH = LH + D2z(lp1+2*lp3, X, Y, Z) +V/C;
    
    M = C*LH;   
    n = 1;  
    
    
function [M, n] = lutt3D_111(inp, X, Y, Z)
    disp('Build 3D Lutt. Hamiltonian (growth:Z[1 1 1] X[1 1 -2] Y[-1 1 0]) (approximate single & mixed derivatives)')
    disp('G. Fishman PRB 52 11132 (1995)')
    %070216 -using ki = -(idxi)   [kij already symmetrized]
    
    conc = inp.lrpot.C;

    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    mparam = inp.mparam; physconst
    
    [lp1, lp2, lp3] = luttingerParameters(conc, mparam, 'INV');
    
    if isfield(inp.VB, 'spherical_approx')
        if inp.VB.spherical_approx  lp3 = (lp2+lp3)/2; lp2 = lp3; disp('----> Spherical band approximaion'); end
    end

    V = bandAlignment(conc, mparam, 'VB')*eV;
    
    if(isfield(inp.VB, 'Epot'))
       if inp.VB.Epot disp('Electric field found: Add electric potential'); end
       V = V + inp.VB.Epot*eV;
    end
    
    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
    C = hbar^2/(me*(1E-9)^2);

    HH =      (1/2) * k2x(lp1+lp3, X, Y, Z);
    HH = HH + (1/2) * k2y(lp1+lp3, X, Y, Z);
    HH = HH + (1/2) * k2z(lp1-2*lp3, X, Y, Z) +V/C;
    
    LH =      (1/2) * k2x(lp1-lp3, X, Y, Z);
    LH = LH + (1/2) * k2y(lp1-lp3, X, Y, Z);
    LH = LH + (1/2) * k2z(lp1+2*lp3, X, Y, Z) +V/C;
    
    R =    -sqrt(3)/2*(          k2x( (2*lp3+lp2)/3 , X, Y, Z) -        k2y( (2*lp3+lp2)/3 , X, Y, Z) - (2i)*kxy( (2*lp3+lp2)/3 , X, Y, Z)   );
    R = R  -4/sqrt(6)*(          kxz(   (lp3-lp2)/2 , X, Y, Z) + (1i) * kyz(   (lp3-lp2)/2 , X, Y, Z)   );
    
    S =       sqrt(3)*(          kxz( (lp3+2*lp2)/3 , X, Y, Z) - (1i) * kyz( (lp3+2*lp2)/3 , X, Y, Z)   );
    S = S  +sqrt(2/3)*(          k2x(   (lp3-lp2)/2 , X, Y, Z) -        k2y(   (lp3-lp2)/2 , X, Y, Z) + (2i)*kxy(   (lp3-lp2)/2 , X, Y, Z)   );
    
    O = sparse(length(HH), length(HH));
    
    if isfield(inp.VB, 'coupling_off')
        if inp.VB.coupling_off  S = O; R = O; disp('----> Valence band coupling switched off'); end
    end
    
    M = C*[HH  -S  R   O;
           -S' LH  O   R;
            R'  O  LH  S;
            O   R' S' HH];   
    n = 4;

      
    
function [M, n] = lutt2D_kz_111(inp, X, Y, Z)
    disp('Build 2D Lutt. Hamiltonian (growth:Z[111]:kz X[11-2] Y[1-10]) (approximate single & mixed derivatives)')
    disp('U. Ekenberg Semicond. Sci. Technol. 4 904-909 (1989) / G. Fishman PRB 52 11132 (1995)')
    %060518
    
    conc = inp.lrpot.C;
 
    [Sy, Sx, Sz] = size(X);
    Stot = Sx*Sy*Sz;
    mparam = inp.mparam; physconst
    
    [lp1, lp2, lp3] = luttingerParameters(conc, mparam, 'INV');
    LP1 = spdiags(reshape(lp1, Stot, 1), 0, Stot, Stot);
    LP2 = spdiags(reshape(lp2, Stot, 1), 0, Stot, Stot);
    LP3 = spdiags(reshape(lp3, Stot, 1), 0, Stot, Stot);
    
    if isfield(inp.VB, 'spherical_approx')
        if inp.VB.spherical_approx  lp3 = (lp2+lp3)/2; lp2 = lp3; disp('----> Spherical band approximaion'); end
    end

    %d/dz = i*kz
    if isfield(inp.VB, 'kz') ikz = 1i*inp.VB.kz; disp(['kz = ' num2str(ikz/(1i)) ' nm-1']); else ikz = 0; end
        
    
    V = bandAlignment(conc, mparam, 'VB')*eV;
    V = spdiags(reshape(V, Stot, 1), 0, Stot, Stot);
    C = -hbar^2/(2*me*(1E-9)^2);

    HH =      D2x(lp1+lp3, X, Y, Z);
    HH = HH + D2y(lp1+lp3, X, Y, Z);
    HH = HH + ikz*ikz*(LP1-2*LP3) + V/C;
        
    LH =      D2x(lp1-lp3, X, Y, Z);
    LH = LH + D2y(lp1-lp3, X, Y, Z);
    LH = LH + ikz*ikz*(LP1+2*LP3) + V/C;
    
    B =   + 2i/sqrt(3)*(       (Dx(X, Y, Z)*(2*LP2+LP3)*ikz + ikz*(2*LP2+LP3)*Dx(X, Y, Z))/2);
    B = B + 2i/sqrt(3)*( (-1i)*(Dy(X, Y, Z)*(2*LP2+LP3)*ikz + ikz*(2*LP2+LP3)*Dy(X, Y, Z))/2);
    B = B - 2i/sqrt(6)*(        D2x(lp2-lp3, X, Y, Z)       - D2y(lp2-lp3, X, Y, Z)) ;
    B = B - 2i/sqrt(6)*( (+1i)*(Dxy(lp2-lp3, X, Y, Z)       + Dyx(lp2-lp3, X, Y, Z))) ;
        
    R =   - 1/sqrt(3)*(          D2x(lp2+2*lp3, X, Y, Z)     - D2y(lp2+2*lp3, X, Y, Z));
    R = R - 1/sqrt(3)*(   (-1i)*(Dxy(lp2+2*lp3, X, Y, Z)     + Dyx(lp2+2*lp3, X, Y, Z)));
    R = R + 4/sqrt(6)*(         (Dx(X, Y, Z)*(LP2-LP3)*ikz   + ikz*(LP2-LP3)*Dx(X, Y, Z))/2);
    R = R + 4/sqrt(6)*(   (+1i)*(Dy(X, Y, Z)*(LP2-LP3)*ikz   + ikz*(LP2-LP3)*Dy(X, Y, Z))/2);
    
    %B = B / (-1i); disp('Transform to U. Ekenberg -> G. Fishman');
    
    O = sparse(length(HH), length(HH));
    
    if isfield(inp.VB, 'coupling_off')
        if inp.VB.coupling_off  B = O; R = O; disp('----> Valence band coupling switched off'); end
    end
    
    M = C*[HH   B  R   O;
            B' LH  O   R;
            R'  O  LH -B;
            O   R' -B' HH];   
    n = 4;
    
function TEST_MODULE(X, Y, Z)
    %%%%%%%%%TEST MODULE: TEST ARBITRARY SET OF DIFF OPERATORS ON ARBITRARY FUNCTION F(X, Y, Z)
    R = sqrt(X.^2+Y.^2+Z.^2);
    %Arbitrary funtion (->0 @ edges)
    F0 = -sin(Z).*cos(X/3).*exp(-R.^2/100);
    F0_sp = spdiags(reshape(F0, Sx*Sy*Sz, 1), 0, Sx*Sy*Sz, Sx*Sy*Sz);

    %Apply arbitrary differential operators
    F = pi*diffy(F0, 2)/(1^2)+diffz(F0, 2)/(1^2) + diffy(F0.*diffz(F0, 1))/1;

    %Define corresponding matrix to test
    M = pi*D2y(ones(size(F0)), X, Y, Z)+D2z(ones(size(F0)), X, Y, Z)+Dyapprx(X, Y, Z)*F0_sp*Dzapprx(X, Y, Z);

    %Solve equation with operator
    F1 = reshape(F, Sx*Sy*Sz, 1);
    F2 = M\F1;
    Fout = reshape(F2, Sy, Sx, Sz);

    %Compare solution (F2) with the input funtion (F)

    %figure; plot(F1, 'k'); hold on; plot(F2, 'b');
    %figure; imagesc([F(:,:, 10) F0(:,:, 10) Fout(:,:, 10)]);daspect([1 1 1]);
    figure; imagesc([reshape(F0(:,1, :), 41, 41) reshape(Fout(:,1, :), 41, 41) reshape(F0(:,1, :), 41, 41)-reshape(Fout(:,1, :), 41, 41)]);daspect([1 1 1]);
%%%%%%%%%END OF TEST MODULE
