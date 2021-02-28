

tic

%Nex=1;

  %  fullenergyN = fullenergy(Q,Nex);
    
    %fullenergyN

%Xxbonding=zeros(Nex);

for Nex=1:2
    
for k=1:Nex
    for h=1:Nex

       
        fullenergyN = fullenergy(Q,Nex);
        
        fullenergynokh = fullenergy(Q,Nex,k,h);
        
        khbonding(k,h,Nex) = fullenergyN - fullenergynokh;

    end
end

end
toc

function Energy = fullenergy(Q,Nex,k,h)

if nargin<3
    k = 0;
    h = -Nex;
end

%disp('!! Check coorinates for Coulomb')
[HH, LH] = comp_luttWF(Q.VB.WF);

WFE = Q.CB.WF;
WFH = sqrt(abs(HH).^2+abs(LH).^2);

%WFH =Q.VB.WF;
WFEN = WFE(:,:,:,1:Nex);

WFHN=WFH(:,:,:,1:Nex);

WF=cat(4,WFEN,WFHN);

h = h+size(WFEN,4);
%size(WF, 4);
%WF=zeros(size(WFE,1),size(WFE,2),size(WFE,3),2*size(WFE,4));
 Eq=zeros(size(WF, 4));

%VBnorm=sqrt(sum(sum(sum(sum(abs(WFH).^2)))));
%CBnorm=sqrt(sum(sum(sum(sum(abs(Q.CB.WF).^2)))));
%WFH=WFH/VBnorm;
%WFE=abs(Q.CB.WF/CBnorm);

for i=1:size(WF, 4);
    
    if (i~=k)&&(i~=h)
    
    for j=1:size(WF, 4);
        
        
        if (j~=i)&&(j~=k)&&(j~=h)
                
        %disp([i j])
    Eq(i, j) = coulomb(WF(:,:,:,i), WF(:,:,:,j), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
    
    
        %Ex(i, j) = N*coulomb(WFE(:,:,:,i), WFH(:,:,:,j), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
        %Ex(i, j) = Occupation(i,j,[Q.CB.E],[Q.VB.E])*coulomb(WFE(:,:,:,i), WFH(:,:,:,j), Q.lrpot.x, Q.lrpot.y, Q.lrpot.z);
        else
            
        end
        
    
    end
    
    else
    
    end
end


Energy=0.5*sum(sum(Eq(:,:)));

end

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

