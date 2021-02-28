function [C, QD1, QD2] = coupled_cylQDs(x, y, z, shape)


[X, Y, Z] = meshgrid(x, y, z);
C = zeros(size(X)); QD1 = C; QD2 = C;

if isfield(shape, 'vqwr_r1')
    R = sqrt(X.^2+Y.^2);
    C(find((R<=shape.vqwr_r1)&(Z>shape.z1))) = 1;
end

if isfield(shape, 'vqwr_r2')
    R = sqrt(X.^2+Y.^2);
    C(find((R<=shape.vqwr_r2)&(Z<shape.z1)&(Z>shape.z2))) = 1;
end

if isfield(shape, 'vqwr_r3')
    R = sqrt(X.^2+Y.^2);
    C(find((R<=shape.vqwr_r3)&(Z<shape.z2))) = 1;
end

if isfield(shape, 'r1')
    R = sqrt(X.^2+Y.^2);
    H = abs(Z-shape.z1);
    C(find((R<shape.r1)&(H<(shape.h1/2)))) = 1;
    QD1(find((R<shape.r1)&(H<(shape.h1/2)))) = 1;

end

if isfield(shape, 'r2')
    R = sqrt(X.^2+Y.^2);
    H = abs(Z-shape.z2);
    C(find((R<shape.r2)&(H<(shape.h2/2)))) = 1;
    QD2(find((R<shape.r2)&(H<(shape.h2/2)))) = 1;
    
end
