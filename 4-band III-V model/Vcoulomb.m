

for k=1:6
    for h=1:6

       
       
        
        VcoulombEx3(k,h) = Vcoulomb2(Q,k,h);
        
       
        

    end
end




function Energy = Vcoulomb2(Q,k,h)



%disp('!! Check coorinates for Coulomb')
[HH, LH] = comp_luttWF(Q.VB.WF);

WFE = Q.CB.WF;
WFH = sqrt(abs(HH).^2+abs(LH).^2);


        
       
       
    Energy = coulomb2(WFE(:,:,:,k), WFH(:,:,:,h), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
    
           

%biEx=Eq(1, 1)+Eq(2, 2)-3*Eq(1, 3);

end


function I = coulomb2(WF1, WF2, x, y, z)

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



