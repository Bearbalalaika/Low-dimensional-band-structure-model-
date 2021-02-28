function [C, QD1, QD2] = ccoupled_ballQDs(x, y, z, shape)

rel1 = 1;
rel2 = 1;
rel3 = 1;

if isfield(shape, 'rel1')
    rel1 = shape.rel1;
end
if isfield(shape, 'rel2')
    rel2 = shape.rel2;
end
if isfield(shape, 'rel3')
    rel3 = shape.rel3;
end

[X, Y, Z] = meshgrid(x, y, z);
C = zeros(size(X)); QD1 = C; QD2 = C;

if isfield(shape, 'vqwr_r1')
    R = sqrt(X.^2+Y.^2);
    C(find((R<=shape.vqwr_r1)&(Z>shape.z1))) = rel1;
end

if isfield(shape, 'vqwr_r2')
    R = sqrt(X.^2+Y.^2);
    C(find((R<=shape.vqwr_r2)&(Z<shape.z1)&(Z>shape.z2))) = rel2;
end

if isfield(shape, 'vqwr_r3')
    R = sqrt(X.^2+Y.^2);
    C(find((R<=shape.vqwr_r3)&(Z<shape.z2))) = rel3;
end

if isfield(shape, 'r1')
    R = sqrt(X.^2+Y.^2+(Z-shape.z1).^2);
    C(find(R<shape.r1)) = 1;
    QD1(find(R<shape.r1)) = 1;
end

if isfield(shape, 'r2')
    R = sqrt(X.^2+Y.^2+(Z-shape.z2).^2);
    C(find(R<shape.r2)) = 1;
    QD2(find(R<shape.r2)) = 1;
end
