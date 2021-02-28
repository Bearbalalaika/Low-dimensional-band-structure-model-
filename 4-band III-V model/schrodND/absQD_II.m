%WFE = 1;
%hh = 0; lh = 1; WFH = []; WFH(1,1,1,1) = hh; WFH(1,1,2, 1) = lh; WFH(1,1,3, 1) = -lh; WFH(1,1,4, 1) = hh; WFH(1,1,1,2) = hh; WFH(1,1,2, 2) = lh; WFH(1,1,3, 2) = lh; WFH(1,1,4, 2) = -hh;


%Verified and correct 060501 !!      See PRB 45 1688 (1992)
function [m, E] = absQD(WFE, WFH, e, EE, EH)

m = [];
if length(e) == 2
    fi = e(1);
    v = e(2);
    BASTARD = 1;
else   
    e = e/sqrt(e(1)^2+e(2)^2+e(3)^2);
    BASTARD = 0;
end

VBnorm=sqrt(sum(sum(sum(sum(abs(WFH).^2)))));
CBnorm=sqrt(sum(sum(sum(sum(abs(WFE).^2)))));
WFH=WFH/VBnorm;
WFE=WFE/CBnorm;

for i=1:size(WFE, 4);
    for j=1:size(WFH, 4);
        if BASTARD
            m(i, j) = M2(WFE(:,:,:,i), WFH(:,:,:,j), fi, v);
        else
            
 % occupation unmber taken from experimental data   
 Pow=40; % nominal excitation power [uWat]
%T=90*0.0862; % exciton effective temperature [meV]
%u=1.561; % effective fermi energy [eV]
deltaQ=41; % differentce of Gs energy position [meV] experiment-theory
Ebg=1.519; % band gap energy [eV]
u=(7.369e-07)*Pow^2+(-3.732e-06)*Pow+(1.549);
T=((0.00144)*Pow^2+(0.4732)*Pow+(5.188))*0.0862;
ucb=2*(u-1.551)/3;
uvb=1*(u-1.551)/3;
Ncb=1/(exp(((EE(i)-EE(1))/1000-ucb)/(T*0.001))+1);
Nvb=1/(exp(((EH(j)-EH(1))/1000-uvb)/(T*0.001))+1);

m(i, j) = Ncb*Nvb*nM2(WFE(:,:,:,i), WFH(:,:,:,j), e);
%N=1/(exp(((EE(i)+EH(j)-deltaQ)/1000+Ebg-u)/(T*0.001))+1);
            
           % m(i, j) = N*nM2(WFE(:,:,:,i), WFH(:,:,:,j), e);
        end
        E(i, j) = EE(i) + EH(j);
    end
end

function m = nM2(E, H, e)

%Basis: [Xu Yu Zu ; Xd Yd Zd]  (u = spin up, d = spin down)

%Valence band basis
VB(:,:,1) = +[1 +i 0; 0 0 0]/sqrt(2);  %3/2
VB(:,:,2) = +[0 0 0; 1 +i 0]/sqrt(6) - [0 0 1; 0 0 0]*sqrt((2/3)); %+1/2
VB(:,:,3) = -[1 -i 0; 0 0 0]/sqrt(6) - [0 0 0; 0 0 1]*sqrt((2/3)); %-1/2
VB(:,:,4) = -[0 0 0; 1 -i 0]/sqrt(2);  %-3/2

%Conduction band
CB(:,:,1) = i*kron([1; 0], [1 1 1]);
CB(:,:,2) = i*kron([0; 1], [1 1 1]);

sz = size(E, 3);
%E = E / sqrt(sum(sum(sum(abs(E(:,:,:)).^2))));
%H = H / sqrt(sum(sum(sum(abs(H(:,:,:)).^2))));

%Integrals
I(1) = sum(sum(sum(E.*H(:,:,(1:sz)))));        %I+3/2
I(2) = sum(sum(sum(E.*H(:,:,sz+(1:sz)))));     %I+1/2
I(3) = sum(sum(sum(E.*H(:,:,2*sz+(1:sz)))));   %I-1/2
I(4) = sum(sum(sum(E.*H(:,:,3*sz+(1:sz)))));   %I-3/2

m = 0;
for cb = 1:2
    m00 = 0;
        for vb=1:4
                m00 = m00 + I(vb)*sum(e.*sum(CB(:,:,cb).*VB(:,:,vb)));
        end
    m = m + abs(m00)^2;
end    



%function m = M2(E, H, fi, v)
%sz = size(E, 3);
%H = H / sqrt(sum(sum(sum(abs(H(:,:,:)).^2))));
%IA = sum(sum(sum(E.*H(:,:,(1:sz)))));        %I+3/2
%IB = sum(sum(sum(E.*H(:,:,sz+(1:sz)))));     %I+1/2
%IC = sum(sum(sum(E.*H(:,:,2*sz+(1:sz)))));   %I-1/2
%ID = sum(sum(sum(E.*H(:,:,3*sz+(1:sz)))));   %I-3/2

%if (v == 0)
%    m = 2/3*(abs(IB).^2+abs(IC).^2);
%elseif (v == pi/2)
%    m = 1/2*(abs(IA).^2+abs(ID).^2) + 1/6*(abs(IB).^2+abs(IC).^2)-cos(2*fi)/sqrt(3) * (IA*conj(IC) + IB*conj(ID));% - 1/sqrt(3) * (IA*IC + IB*ID) * cos(2*pi);
%else
%    m = NaN;
%end

