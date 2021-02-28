function P = pyralens;
x = -6:6
C = ones(size(x));
C(1) = 0; C(length(x)) = 0;
P = zeros(size(C));
P(2) = 1;
P(11) = 2;
figure;
for i = 1:1000
%P = (shift(P, -1) + shift(P, 1) + shift(P, -2) + shift(P, 2))/5;
P = (shift(P, -2) + shift(P, 2) )/2;
P(1) = P(2)-1;
P(13) = P(12)-2;
%if ~mod(i, 10) plot(P);  f = getframe(gcf); disp(i); end;
P = P-max(P);
end

figure; plot(x, P)


function A = shift(A, dim)
[SY, SX, SZ] = size(A);
for i=1:length(dim)
if (dim(i) == -1)
    A(2:SY,:,:) = A(1:(SY-1),:,:);
elseif (dim(i) == 1)
    A(1:(SY-1),:,:) = A(2:SY,:,:);
elseif (dim(i) == -2)
    A(:,2:SX,:) = A(:,1:(SX-1),:);
elseif (dim(i) == 2)
    A(:,1:(SX-1),:) = A(:,2:SX,:); 
elseif (dim(i) == -3)
    A(:,:,2:SZ) = A(:,:,1:(SZ-1));
elseif (dim(i) == 3)
    A(:,:,1:(SZ-1)) = A(:,:,2:SZ);
end
end
