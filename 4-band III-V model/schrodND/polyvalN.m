function v=polyvalN(p, x)
lp = length(p); x0 = ones(size(x)); v = zeros(size(x));
for i=1:lp
    v = v + x0*p(lp-i+1);
    x0 = x0.*x;
end
