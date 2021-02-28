%[Fg, G] = gaussBlur(F, dx, dy, dz)
%F = funtion to filter
%dx, dy, dz = Gaussian FWHM
%Fg = Filtered function
%G = Convolution kernel

function [Fg, G] = gaussBlur(F, dx, dy, dz)

x = -1:1; y = x; z = x;
w = 3;
disp('Gaussian convolution')

if (dx(1)~=0)
    if dx>=1; x = -ceil(w*dx):ceil(w*dx); end
    Gx = analvoigt(1, 0, 0, dx(1), x); Gx = Gx/sum(Gx);
    Fg = gaussBlur_x(F, Gx);
else 
    Fg = F;
end


if nargin>2
    if (dy(1)~=0)
        if dy>=1; y = -ceil(w*dy):ceil(w*dy); end
        Gy = analvoigt(1, 0, 0, dy(1), y); Gy = Gy/sum(Gy);
        Fg = gaussBlur_y(Fg, Gy);
    else 
        Fg = Fg;
    end
end

if nargin>3
    if (dz(1)~=0)
        if dz>=1; z = -ceil(w*dz):ceil(w*dz); end
        Gz = analvoigt(1, 0, 0, dz(1), z); Gz = Gz/sum(Gz);
        Fg = gaussBlur_z(Fg, Gz);
    else 
        Fg = Fg;
    end
end

function Fg = gaussBlur_x(F0, G)
    [sy, sx, sz] = size(F0);
    N = (length(G)-1)/2;
    F = zeros(sy, sx+2*N, sz);
    for i=1:N
        F(:, i, :) = F0(:,1,:);
        F(:, N+sx+(i), :) = F0(:,sx,:);
    end

    for i=1:sx;
        F(:, N + (i), :) = F0(:,i,:);
    end
    
    G = reshape(G, 1, length(G), 1);
    Fg = convn(F, G, 'valid');

function Fg = gaussBlur_y(F0, G)
    [sy, sx, sz] = size(F0);
    N = (length(G)-1)/2;
    F = zeros(sy+2*N, sx, sz);
    for i=1:N
        F(i, :, :) = F0(1,:,:);
        F(N+sy+(i), :, :) = F0(sy,:,:);
    end
    
    for i=1:sy
        F(N + (i), :, :) = F0(i,:,:);
    end
    
    G = reshape(G, length(G),1, 1);
    Fg = convn(F, G, 'valid');

function Fg = gaussBlur_z(F0, G)
    [sy, sx, sz] = size(F0);
    N = (length(G)-1)/2;
    F = zeros(sy, sx, sz+2*N);
    for i=1:N
        F(:, :, i) = F0(:,:,1);
        F(:, :, N+sz+(i)) = F0(:,:,sz);
    end
    
    for i=1:sz
        F(:, :, N + (i)) = F0(:,:,i);
    end
    
    G = reshape(G, 1, 1, length(G));
    Fg = convn(F, G, 'valid');
