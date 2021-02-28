%C = resampleND(xh, yh, zh, C, xl, yl, zl)
%Resample C defined on 1D coordinates xh, yh, zh
%on new grid defined by xl, yl, zl;
%(zl can exceed zh, zl>zh = max(zh))

function C = resampleND(xh, yh, zh, C, xl, yl, zl)

[Sy, Sx, Sz] = size(C);

if (Sx>1) dxh = mean(diff(xh)); else dxh = 1; end
if (Sy>1) dyh = mean(diff(yh)); else dyh = 1; end
if (Sz>1) dzh = mean(diff(zh)); else dzh = 1; end

if (Sx>1) dxl = mean(diff(xl)); else dxl = 0; end
if (Sy>1) dyl = mean(diff(yl)); else dyl = 0; end
if (Sz>1) dzl = mean(diff(zl)); else dzl = 0; end

[Xh, Yh, Zh] = meshgrid(xh, yh, zh);
[Xl, Yl, Zl] = meshgrid(xl, yl, zl);
%Zl(find(Zl>max(zh))) = max(zh);
%Zl(find(Zl<min(zh))) = min(zh);

C = gaussBlur(C, 0.85*dxl/dxh, 0.85*dyl/dyh, 0.85*dzl/dzh);


dim = ~~(Sx-1) + ~~(Sy-1) + ~~(Sz-1);

if (dim == 3)
    C = interp3(Xh, Yh, Zh, C, Xl, Yl, Zl);
end

if (dim == 2)
    if Sy>1
        C = interp2(Xh, Yh, C, Xl, Yl);
    else
        C0 = interp2(squeeze(Xh)', squeeze(Zh)', squeeze(C)', squeeze(Xl)', squeeze(Zl)');
        [Sy, Sx, Sz] = size(Xl);
        C = zeros(1, Sx, Sz);
        for i=1:Sz
            C(1,:,i) = C0(i,:);
        end
    end
end
