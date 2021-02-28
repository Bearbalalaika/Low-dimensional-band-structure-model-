function C = shiftZ(C, n)
sz = size(C, 3);
if n>0
    C0 = C(:,:,1);
else
    C0 = C(:,:,sz);
end

C = circshift(C, [0 0 1]*n);

for i=1:n
    if n>0
        C(:,:,i) = C0;
    else
        C(:,:,sz-(i+1)) = C0;
    end
end
