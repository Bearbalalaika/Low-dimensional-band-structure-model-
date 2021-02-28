%Rot 90 deg x-y also for 3D arrays
function C = rot90ND(C);
[sy, sx, sz] = size(C);
for i=1:sz;
    C(:,:,i) = rot90(C(:,:,i));
end