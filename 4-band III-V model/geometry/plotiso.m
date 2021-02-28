function plotiso(C, WF, v, color)
WF = abs(WF).^2;

if nargin<4
    color = [0 0 1]
end

if length(C)>1
p = patch(isosurface(C, 0.1)); isonormals(C, p);
set(p, 'FaceColor', 0.5*[1 1 1], 'EdgeColor', 'none');
alpha(p, 0.33);
end

p = patch(isosurface(WF, v*0.25)); isonormals(WF, p);
set(p, 'FaceColor', color, 'EdgeColor', 'none');
alpha(p, 0.6);


%p = patch(isosurface(WF, v)); isonormals(WF, p);
%set(p, 'FaceColor', color/2, 'EdgeColor', 'none');
%alpha(p, 0.6);


daspect([1 1 1]); axis tight;

%camup([1 0 0 ]); campos([25 -55 5]) 
camlight; lighting gouraud