function Q = exciton_pert(Q)
disp('!! Check coorinates for Coulomb')
[HH, LH] = comp_luttWF(Q.VB.WF);
WFE = Q.CB.WF;
WFH = sqrt(abs(HH).^2+abs(LH).^2);
%WFH =Q.VB.WF;

%VBnorm=sqrt(sum(sum(sum(sum(abs(WFH).^2)))));
%CBnorm=sqrt(sum(sum(sum(sum(abs(Q.CB.WF).^2)))));
%WFH=WFH/VBnorm;
%WFE=abs(Q.CB.WF/CBnorm);


for i=1:size(WFE, 4);
    for j=1:size(WFH, 4);
        disp([i j])
%{    
 % occupation unmber taken from experimental data     
%T=36*0.0862; % exciton effective temperature [meV]
%u=1.551; % effective fermi energy [eV]
deltaQ=41; % differentce of Gs energy position [meV] experiment-theory
Ebg=1.519; % band gap energy [eV]
%N=1/(exp(((Q.CB.E(i)+Q.VB.E(j)-deltaQ)/1000+Ebg-u)/(T*0.001))+1);
         Pow=20; % nominal excitation power [uWat]
%T=90*0.0862; % exciton effective temperature [meV]
%u=1.561; % effective fermi energy [eV]
deltaQ=41; % differentce of Gs energy position [meV] experiment-theory
Ebg=1.519; % band gap energy [eV]
u=(7.369e-07)*Pow^2+(-3.732e-06)*Pow+(1.549);
T=((0.00144)*Pow^2+(0.4732)*Pow+(5.188))*0.0862;
ucb=2*(u-1.551)/3;
uvb=(u-1.551)/3;
Ncb=1/(exp(((Q.CB.E(i)-Q.CB.E(1))/1000-ucb)/(T*0.001))+1);
Nvb=1/(exp(((Q.VB.E(j)-Q.VB.E(1))/1000-uvb)/(T*0.001))+1);
%}
Ex(i, j) = coulomb(WFE(:,:,:,i), WFH(:,:,:,j), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
        %Ex(i, j) = N*coulomb(WFE(:,:,:,i), WFH(:,:,:,j), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
        %Ex(i, j) = Occupation(i,j,[Q.CB.E],[Q.VB.E])*coulomb(WFE(:,:,:,i), WFH(:,:,:,j), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
    end
end
Q.X.info = 'Perturbation';
Q.X.Ex = kron(Ex, [1 1]);



function I = coulomb(WF1, WF2, x, y, z)

%Coulomb coordinatesQ
dxyz_c = 1;
xc = -30:dxyz_c:30;
yc = -30:dxyz_c:30;
zc = -150:dxyz_c:150;

WF1 = resampleND(x, y, z, WF1, xc, yc, zc); WF1(isnan(WF1)) = 0; %Conserve density
WF1 = WF1/sqrt(sum(abs(WF1(:)).^2)); %Normalize

WF2 = resampleND(x, y, z, WF2, xc, yc, zc); WF2(isnan(WF2)) = 0; %Conserve density
WF2 = WF2/sqrt(sum(abs(WF2(:)).^2)); %Normalize

I = FFTn_coulomb(dxyz_c, WF1, WF2, 'integral');
