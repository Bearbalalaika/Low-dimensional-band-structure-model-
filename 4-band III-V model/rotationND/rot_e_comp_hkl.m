%Rotate coordinate system for alla strain components
%rot_e_comp(i, hkl, NORMorINV)
%Verified, works perfectly 050417

function o=rot_e_comp(i, hkl, NORMorINV)

if isnan(hkl)
    Q = eye(3, 3);
else
    Q = rotMatrix_hkl(hkl);
end

if strcmp(NORMorINV, 'INV') 
    disp('INVERSE transformation of strain components')
    Q = inv(Q);
end

Q11 = Q(1,1); Q12 = Q(1,2); Q13 = Q(1,3);
Q21 = Q(2,1); Q22 = Q(2,2); Q23 = Q(2,3);
Q31 = Q(3,1); Q32 = Q(3,2); Q33 = Q(3,3);


%for the following to work: NEW_R = Q*OLD_R

xx = i.e_xx;
yy = i.e_yy;
zz = i.e_zz;
yz = i.e_yz/2;
xz = i.e_xz/2;
xy = i.e_xy/2;

%Strain in New coordinates 
nxx = (Q11*xx+Q12*xy+Q13*xz)*(Q11)+(Q11*xy+Q12*yy+Q13*yz)*(Q12)+(Q11*xz+Q12*yz+Q13*zz)*(Q13);
nxy = (Q11*xx+Q12*xy+Q13*xz)*(Q21)+(Q11*xy+Q12*yy+Q13*yz)*(Q22)+(Q11*xz+Q12*yz+Q13*zz)*(Q23);
nxz = (Q11*xx+Q12*xy+Q13*xz)*(Q31)+(Q11*xy+Q12*yy+Q13*yz)*(Q32)+(Q11*xz+Q12*yz+Q13*zz)*(Q33);
nyy = (Q21*xx+Q22*xy+Q23*xz)*(Q21)+(Q21*xy+Q22*yy+Q23*yz)*(Q22)+(Q21*xz+Q22*yz+Q23*zz)*(Q23);
nyz = (Q21*xx+Q22*xy+Q23*xz)*(Q31)+(Q21*xy+Q22*yy+Q23*yz)*(Q32)+(Q21*xz+Q22*yz+Q23*zz)*(Q33);
nzz = (Q31*xx+Q32*xy+Q33*xz)*(Q31)+(Q31*xy+Q32*yy+Q33*yz)*(Q32)+(Q31*xz+Q32*yz+Q33*zz)*(Q33);

%Create output strucure
o = i;
o.comment = ['ROTATED  [' num2str(hkl) '] ' NORMorINV];
o.e_xx = nxx;
o.e_yy = nyy;
o.e_zz = nzz;
o.e_yz = nyz*2;
o.e_xz = nxz*2;
o.e_xy = nxy*2;
