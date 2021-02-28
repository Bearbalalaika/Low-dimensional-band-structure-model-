%Mout = resize(Min, factor [, method]);
%
%Min = matrix in
%factor = resize facor (0.5 = half)
%method = interpolation method (default linear)

function [Mout, sizeOut] = resize(Min, factor, method);
factor
if (length(factor) == 1) sizeOut = round(size(Min)*factor)
else sizeOut = factor;
end
Mout = resizeND(Min, sizeOut(2), sizeOut(1));

disp('No method implemented')

%Old version...
%S = size(Min);
%
%if length(S) == 3
%   for i=1:S(3)
%      [Mout(:,:,i), sizeOut] = resize(Min(:,:,i), factor, method);
%   end
%else
%
%if ((min(factor) == 1)& (max(factor) == 1))
%   Mout = Min;
%   sizeOut = size(Mout);
%else
%   
%
%sy = S(1);
%sx = S(2);

%if length(factor)==1
%   sizeOut = round(factor * [sy sx]);
%else
%   sizeOut = factor;
%end
%
%   if nargin < 3
%   	method = 'linear';
%   end
%   
%   disp('Filtering disabled in resize()');
%   if min(factor) < 0.75
%      Min = filter2([1 2 1; 2 4 2; 1 2 1]/16,Min);
%   end
%   
%	Min = double(Min);
%
%	x_out = 1:sizeOut(2);
%   y_out = 1:sizeOut(1);
   
%   x_int = 1 + (x_out-1) * (sx-1)/x_out(length(x_out)-1);
%   y_int = 1 + (y_out-1) * (sy-1)/y_out(length(y_out)-1);
%   oneVector = ones(size(x_int));
%   for i=1:length(y_out)
%      X(i,:) = x_int;
%      Y(i,:) = y_int(i)*oneVector;
%   end
   
%   Mout = interp2(Min, X, Y, method);
%end
%end
