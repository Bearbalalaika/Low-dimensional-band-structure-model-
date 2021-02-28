function A = mirror(A);
[sy, sx] = size(A); B(:,1:sx) = A(:,abs(-sx:-1));
A = B;
