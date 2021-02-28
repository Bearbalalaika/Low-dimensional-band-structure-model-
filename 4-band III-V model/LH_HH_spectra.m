function [ ml, mh, E ] = LH_HH_spectra( WFE, WFH, EE, EH )
%LH_HH_SPECTRA Summary of this function goes here

% This script plots matrix elements of dipol matrix which correspond to HH
% and LH tansitions


ml = [];
mh = [];
 


VBnorm=sqrt(sum(sum(sum(sum(abs(WFH).^2)))));
CBnorm=sqrt(sum(sum(sum(sum(abs(WFE).^2)))));
WFH=WFH/VBnorm;
WFE=WFE/CBnorm;

for k=1:size(WFE, 4)
    for j=1:size(WFH, 4)
       
%T=36*0.0862; % exciton effective temperature [meV]
%u=1.551; % effective fermi energy [eV]
deltaQ=41; % differentce of Gs energy position [meV] experiment-theory
Ebg=1.519; % band gap energy [eV]
%N=1/(exp(((Q.CB.E(i)+Q.VB.E(j)-deltaQ)/1000+Ebg-u)/(T*0.001))+1);
         Pow=60; % nominal excitation power [uWat]
%T=90*0.0862; % exciton effective temperature [meV]
%u=1.561; % effective fermi energy [eV]
deltaQ=41; % differentce of Gs energy position [meV] experiment-theory
Ebg=1.519; % band gap energy [eV]
u=(7.369e-07)*Pow^2+(-3.732e-06)*Pow+(1.549);
T=((0.00144)*Pow^2+(0.4732)*Pow+(5.188))*0.0862;
ucb=2*(u-1.551)/3;
uvb=(u-1.551)/3;
Ncb=1/(exp(((EE(k)-EE(1))/1000-ucb)/(T*0.001))+1);
Nvb=1/(exp(((EH(j)-EH(1))/1000-uvb)/(T*0.001))+1);

%Ex(i, j) = Ncb*Nvb*coulomb(WFE(:,:,:,i), WFH(:,:,:,j), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
%N=1/(exp(((EE(k)+EH(j)-deltaQ)/1000+Ebg-u)/(T*0.001))+1);
        
            ml(k, j) = Ncb*Nvb*nM2l(WFE(:,:,:,k), WFH(:,:,:,j));
            mh(k, j) = Ncb*Nvb*nM2h(WFE(:,:,:,k), WFH(:,:,:,j));
            
        E(k, j) = EE(k) + EH(j);
    end

end


function m = nM2l(E, H)


sz = size(E, 3);
%E = E / sqrt(sum(sum(sum(abs(E(:,:,:)).^2))));
%H = H / sqrt(sum(sum(sum(abs(H(:,:,:)).^2))));

%Integrals
I(1) = sum(sum(sum(E.*H(:,:,(1:sz)))));        %I+3/2
I(2) = sum(sum(sum(E.*H(:,:,sz+(1:sz)))));     %I+1/2
I(3) = sum(sum(sum(E.*H(:,:,2*sz+(1:sz)))));   %I-1/2
I(4) = sum(sum(sum(E.*H(:,:,3*sz+(1:sz)))));   %I-3/2


m =abs(I(2))^2+abs(I(3))^2;
%m =abs(I(2)+I(3))^2;
end    


function m = nM2h(E, H)



sz = size(E, 3);
%E = E / sqrt(sum(sum(sum(abs(E(:,:,:)).^2))));
%H = H / sqrt(sum(sum(sum(abs(H(:,:,:)).^2))));

%Integrals
I(1) = sum(sum(sum(E.*H(:,:,(1:sz)))));        %I+3/2
I(2) = sum(sum(sum(E.*H(:,:,sz+(1:sz)))));     %I+1/2
I(3) = sum(sum(sum(E.*H(:,:,2*sz+(1:sz)))));   %I-1/2
I(4) = sum(sum(sum(E.*H(:,:,3*sz+(1:sz)))));   %I-3/2


m = abs(I(1))^2+abs(I(4))^2;
%m = abs(I(1)+I(4))^2;


end  

end