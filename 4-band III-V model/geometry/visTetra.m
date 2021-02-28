x0 = -20:0.5:20;
y0 = -20:0.5:20;
z0 = -20:0.5:20;
Sx = length(x0); Sy = length(y0); Sz = length(z0);

C1 = tetra(30);
C2 = tetra(26);
[C3, Ccut] = tetra(20);
C = C1+C2 +1 * C3;
C0 = C;
C = C + 0*Ccut;




%C1 = gaussBlur(sum(C1, 3), 10, 10, 1);
%C3 = gaussBlur(sum(C3, 3), 15, 15, 1);
%Ccut = Ccut(:,:,1);
%C1 = C1 + 0*Ccut;
%C3 = C3 + 0*Ccut;

%figure; h1=surfl(C1); hold on; h2=surfl(C3);
%set(h1, 'facecolor', [1 0 0]); alpha(h1, 0.6)
%set(h2, 'facecolor', [0 0 1]); alpha(h2, 0.6)
%V = zeros(Sy, Sx, Sz);
%for i=1:Sx
    %for j = 1:Sy
      %  V(j, i, round(21+(C3(j, i):C1(j, i)))) = 1;
        %    end
    %end
%size(V)
%size(Ccut)
%figure;
%isosurface(0*Ccut+gaussBlur(V, 3, 3, 3));


%break
figure; isosurface(x0, y0, z0, gaussBlur(C, 2, 2, 2));alpha(0.7);

[x, y, z] = meshgrid(x0, y0, z0);
Ccut = 1 + Ccut; Ccut(isnan(Ccut)) = 0;
C0(find(C0<2)) = NaN;
%C0 = C0 + Ccut;
%x(~isnan(C)) = NaN;
%y(~isnan(C)) = NaN;
%z(~isnan(C)) = NaN;

xs = reshape(x(41,:,:), 81, 81)';
ys = reshape(y(41,:,:), 81, 81)';
zs = reshape(z(41,:,:), 81, 81)';
S = reshape(C0(41,:,:), 81, 81)';
hold on;
s = surf(xs,ys,zs,S);shading flat
alpha(s, 0.5);

[x, y, z] = meshgrid(x0, y0, z0);

%y(~isnan(C)) = NaN;
%z(~isnan(C)) = NaN;

xs = reshape(x(:,42,:), 81, 81)';
ys = reshape(y(:,42,:), 81, 81)';
zs = reshape(z(:,42,:), 81, 81)';
S = reshape(C0(:,42,:), 81, 81)';
hold on;
s = surf(xs,ys,zs,S);shading flat
alpha(s, 0.5);