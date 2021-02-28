%M = finitDiff1D(D2, D1, D0, C2, C1, C0, dz)

function M = finitDiff1D(D2, D1, D0, C2, C1, C0, dz)

SZ = max([length(D2) length(D1) length(D0)]);

if length(C1) == 1
    C1 = C1 * ones(SZ,1);
end

if length(C2) == 1
    C2 = C2 * ones(SZ,1);
end


C1 = C1*[1 1];
C2 = C2*[1 1 1];


i = 1:SZ; i = i';
M0 = spdiags(C0.*D0, 0, SZ, SZ);

M1 = 0.5*spdiags(C1.*[-D1 D1], [-1 1], SZ, SZ);

i = (1:(SZ-1))';
D2z(i,1) = (D2(i) + D2(i+1))/2;
j = (2:(SZ-1))';
M2 = spdiags(C2.*[[D2z; NaN] -[D2(1)+D2z(1); D2z(j-1)+ D2z(j); D2(SZ) + D2z(SZ-1)] [NaN;D2z]], [-1 0 1], SZ, SZ);
M = M2/dz^2 + M1/dz + M0;
