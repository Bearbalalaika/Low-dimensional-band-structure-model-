function [HH, LH, hh, lh] = comp_luttWF(WF)
[Sy, Sx, Sz, n] = size(WF);
LZ = Sz/4;
LH = zeros(Sy, Sx, LZ, n/2);
HH = zeros(Sy, Sx, LZ, n/2);
hh = zeros(n/2, 1);
lh = zeros(n/2, 1);

for i=1:(n/2)
    indx = 1:LZ;
    HH1 = WF(:,:,0*LZ+indx,2*i-1);
    HH2 = WF(:,:,3*LZ+indx,2*i-1);
    HH3 = WF(:,:,0*LZ+indx,2*i);
    HH4 = WF(:,:,3*LZ+indx,2*i);

    LH1 = WF(:,:,1*LZ+indx,2*i-1);
    LH2 = WF(:,:,2*LZ+indx,2*i-1);
    LH3 = WF(:,:,1*LZ+indx,2*i);
    LH4 = WF(:,:,2*LZ+indx,2*i);
    
    
    % /2 ?? Why?
    
    HH(:,:,:,i) = sqrt((abs(HH1).^2 + abs(HH2).^2 + abs(HH3).^2 + abs(HH4).^2)/2);
    LH(:,:,:,i) = sqrt((abs(LH1).^2 + abs(LH2).^2 + abs(LH3).^2 + abs(LH4).^2)/2);
    
    %HH(:,:,:,i) = sqrt((abs(HH1).^2 + abs(HH2).^2 + abs(HH3).^2 + abs(HH4).^2));
    %LH(:,:,:,i) = sqrt((abs(LH1).^2 + abs(LH2).^2 + abs(LH3).^2 + abs(LH4).^2));
    
    
    hh(i) = sum(sum(sum(HH(:,:,:,i).^2)));
    lh(i) = sum(sum(sum(LH(:,:,:,i).^2)));
    
    S = hh(i) + lh(i);

    hh(i) = hh(i)/S;
    lh(i) = lh(i)/S;
    

    
end