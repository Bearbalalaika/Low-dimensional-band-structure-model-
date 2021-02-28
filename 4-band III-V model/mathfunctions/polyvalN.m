function v = polyvalN(p, x);
lp = length(p); v = zeros(size(x)); x0 = ones(size(x));
for i=1:lp
    v = v + p(lp+1-i)*x0;
    x0 = x0.*x;
end
