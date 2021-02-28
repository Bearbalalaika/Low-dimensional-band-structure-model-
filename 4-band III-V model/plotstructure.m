%clear all;
%load Double_QD_2x9nm_barrier_6nm_20_Barrier.mat;
%[HH, LH] = comp_luttWF(Q.VB.WF);
%figure;
%p1=patch(isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, Q.lrpot.C, 0.34));
%isonormals(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, Q.lrpot.C, p1); 
%daspect([1 1 1]);
%set(p1, 'FaceColor', [0 0 0], 'EdgeColor', 'none');
%axis tight;camlight; lighting gouraud;
%alpha(p1,0.2);
%hold on;
figure;
% p2=patch(isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z,Q.lrpot.C, 0.02));
% isonormals(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, Q.lrpot.C, p2); 
% daspect([1 1 1]);
% set(p2, 'FaceColor', [0 0 0], 'EdgeColor', 'none');
% axis tight;camlight; lighting gouraud;
% alpha(p2,0.05);
% hold on;
p3=patch(isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z,Q.lrpot.C, 0.1));
%p3=isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z,Q.lrpot.C, 0.25);

%isonormals(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, Q.lrpot.C, p3); 
daspect([1 1 1]);
set(p3, 'FaceColor', [0 0 0], 'EdgeColor', 'none');



axis tight;
axis([-20 20 -20 20 -30 30]);
%set(gca, 'XTick', [-10; 0; 10],'YTick', [-10; 0; 10],'XTickLabel', [],'YTickLabel', []);
camlight; 
lighting gouraud;
%lighting phong;
%lighting flat;
alpha(p3,0.25);
set(gca, 'FontSize', 26,'Linewidth',3,'TickLength',[0.03 0.025]);
%set(gca, 'XTick', [-10; 10],'YTick', [-10; 10],'ZTick', [-40; -10; 10; 40]);
