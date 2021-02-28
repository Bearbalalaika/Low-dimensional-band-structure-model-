% See PRB 45 1688 (1992)
function [m, E] = absQWR(WFE, WFH, e, EE, EH)

m = [];
if length(e) == 2
    fi = e(1);
    v = e(2);
    BASTARD = 1;
else   
    e = e/sqrt(e(1)^2+e(2)^2+e(3)^2);
    BASTARD = 0;
end


for i=1:size(WFE, 4);
    for j=1:size(WFH, 4);
        if BASTARD
            m(i, j) = M2(WFE(:,:,i), WFH(:,:,j), fi, v);
        else
            m(i, j) = nM2(WFE(:,:,i), WFH(:,:,j), e);
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

sx = size(E, 2);
H = H / sqrt(sum(sum(abs(H(:,:,:)).^2)));

%Integrals
I(1) = sum(sum(E.*H(:,(1:sx))));        %I+3/2
I(2) = sum(sum(E.*H(:,sx+(1:sx))));     %I+1/2
I(3) = sum(sum(E.*H(:,2*sx+(1:sx))));   %I-1/2
I(4) = sum(sum(E.*H(:,3*sx+(1:sx))));   %I-3/2

m = 0;
for cb = 1:2
    m00 = 0;
        for vb=1:4
                m00 = m00 + I(vb)*sum(e.*sum(CB(:,:,cb).*VB(:,:,vb)));
        end
    m = m + abs(m00)^2;
end    



function m = M2(E, H, fi, v)
sx = size(E, 2);
H = H / sqrt(sum(sum(sum(abs(H(:,:,:)).^2))));
IA = sum(sum(sum(E.*H(:,(1:sx)))));        %I+3/2
IB = sum(sum(sum(E.*H(:,sx+(1:sx)))));     %I+1/2
IC = sum(sum(sum(E.*H(:,2*sx+(1:sx)))));   %I-1/2
ID = sum(sum(sum(E.*H(:,3*sx+(1:sx)))));   %I-3/2

if (v == 0)
    m = 2/3*(abs(IB).^2+abs(IC).^2);
elseif (v == pi/2)
    m = 1/2*(abs(IA).^2+abs(ID).^2) + 1/6*(abs(IB).^2+abs(IC).^2)-cos(2*fi)/sqrt(3) * (IA*conj(IC) + IB*conj(ID));% - 1/sqrt(3) * (IA*IC + IB*ID) * cos(2*pi);
else
    m = NaN;
end

