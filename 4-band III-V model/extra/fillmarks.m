function fillmarks(h)
if nargin<1
    h = gco;
end

set(h, 'MarkerFaceColor', get(h, 'color'));
