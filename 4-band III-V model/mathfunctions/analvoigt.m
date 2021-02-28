%value = analvoigt(a_L, pos, gamma_L, gamma_G, x)
%a_L = Lorenzian amplitude
%pos = position
%gamma_L = FWHM Lorenzian
%gamma_G = FWHM Gaussian
%x = coordinate

function V = analvoigt(a_L, pos, gamma_L, gamma_G, x)

if (gamma_L == 0)		%Gaussian
   gamma_G = 0.6*gamma_G;
   V = a_L * exp(-((x-pos)/gamma_G).^2);
elseif (gamma_G == 0)		%Lorenzian
   gamma_L = 0.5*gamma_L;
   V = a_L ./ (1 + ((x-pos)/gamma_L).^2);
elseif ((gamma_G ~= 0) & (gamma_L ~= 0))
   x(length(x)+1) = pos;
	V = zeros(1, length(x));
	sqrtln2 = sqrt(log(2));
	sqrtpi = sqrt(pi);
	X = (x-pos)*2*sqrtln2/gamma_G;
	Y = gamma_L*sqrtln2/gamma_G;
	A = -[1.2150 1.3509 1.2150 1.3509];
	B = [1.2359 0.3786 -1.2359 -0.3786];
	C = [-0.3085 0.5906 -0.3085 0.5906];
	D = [0.0210 -1.1858 -0.021 1.1858];

for i=1:4
     	V = V + (C(i)*(Y-A(i))+D(i)*(X-B(i))) ./ ((Y-A(i))^2 + (X-B(i)).^2);
	end
   V0 = V(length(V));
   V(length(V)) = [];
	V = a_L*V/V0; %*gamma_L*a_L*sqrtpi*sqrtln2/gamma_G;
else   
	V = zeros(1, length(x));   
end
