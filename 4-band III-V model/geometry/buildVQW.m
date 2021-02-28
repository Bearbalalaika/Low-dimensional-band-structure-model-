function C = buildVQW(x, y, z, a, x0, y0)
    [X, Y, Z] = meshgrid(x-x0, y-y0, z);
    L = a;
    C = ones(size(X));

    f = -0*pi/180;
    M = sin(f)*X-cos(f)*Y;
    %sin_f = 1; cos_f = 0;
    %M = sin_f*X-cos_f*Y;
    C(find((X>=0)&(abs(M)<(L/2)))) = 0;

    f = -60*pi/180;
    M = sin(f)*X-cos(f)*Y;
    C(find((Y>=0)&(abs(M)<(L/2)))) = 0;

    f = 60*pi/180;
    M = sin(f)*X-cos(f)*Y;
    C(find((Y<=0)&(abs(M)<(L/2)))) = 0;
    C = 1-C;