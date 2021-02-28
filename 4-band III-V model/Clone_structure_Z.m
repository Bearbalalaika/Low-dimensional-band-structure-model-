function C = Clone_structure_Z(C, cz, sz, Csingel);

if nargin<4
    Csingel = C;
end

C0 =  Csingel(:,:,1);
for i=1:size(C, 3)
    Csingel(:,:,i) = Csingel(:,:,i) - C0;
    C(:,:,i) = C(:,:,i) - C0;
end

for i=1:length(cz)
    C = C + shiftZ(Csingel, cz(i));
end
clear Csingel;

C = shiftZ(C, sz);

for i=1:size(C, 3)
    C(:,:,i) = C(:,:,i) + C0;
end
