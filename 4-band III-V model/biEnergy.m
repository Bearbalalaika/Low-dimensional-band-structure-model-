








%disp('!! Check coorinates for Coulomb')
[HH, LH] = comp_luttWF(Q.VB.WF);

WFE = Q.CB.WF;
WFH = sqrt(abs(HH).^2+abs(LH).^2);

%WFH =Q.VB.WF;




%size(WF, 4);
%WF=zeros(size(WFE,1),size(WFE,2),size(WFE,3),2*size(WFE,4));


%VBnorm=sqrt(sum(sum(sum(sum(abs(WFH).^2)))));
%CBnorm=sqrt(sum(sum(sum(sum(abs(Q.CB.WF).^2)))));
%WFH=WFH/VBnorm;
%WFE=abs(Q.CB.WF/CBnorm);


        
        
   EE1= coulomb(WFE(:,:,:,1), WFE(:,:,:,1), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
                   
   HH1= coulomb(WFH(:,:,:,1), WFH(:,:,:,1), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
   
   
  EQ = coulomb(WFE(:,:,:,1), WFH(:,:,:,1), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
                
         
    
    


%Energy=0.5*sum(sum(Eq(:,:)));


biXEnergy=EE1+HH1-2*EQ

triExE=EE1-EQ

triHxE=HH1-EQ

exiton=EQ

function I = coulomb(WF1, WF2, x, y, z)

%Coulomb coordinatesQ
dxyz_c = 1;
xc = -30:dxyz_c:30;
yc = -30:dxyz_c:30;
zc = -250:dxyz_c:250;

WF1 = resampleND(x, y, z, WF1, xc, yc, zc); WF1(isnan(WF1)) = 0; %Conserve density
WF1 = WF1/sqrt(sum(abs(WF1(:)).^2)); %Normalize

WF2 = resampleND(x, y, z, WF2, xc, yc, zc); WF2(isnan(WF2)) = 0; %Conserve density
WF2 = WF2/sqrt(sum(abs(WF2(:)).^2)); %Normalize

I = FFTn_coulomb(dxyz_c, WF1, WF2, 'integral');
end



