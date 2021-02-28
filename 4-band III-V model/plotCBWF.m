 %clear all;
%load Double_QD_2x9nm_barrier_6nm_20_Barrier.mat;
%hold on;
figure;

p0=patch(isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(Q.CB.WF(:,:,:,1)).^2, 0.0003));
isonormals(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(Q.CB.WF(:,:,:,1)).^2, p0); 
daspect([1 1 1]);
set(p0, 'FaceColor', [0 1 0], 'EdgeColor', 'none');
%axis tight;
%axis ([-15 15 -15 15 -80 80]);
camlight; 
lighting gouraud;
%lighting phong;
%lighting flat;
alpha(p0,0.7);box on;
axis on;
axis ([-10 10 -10 10 -20 20]);
%set(gca, 'XTick', [],'YTick', []);
%set(gca, 'ZTickLabel','','Linewidth',1.5,'TickLength',[0.02,0.02]);
set(gca,'Linewidth',4,'TickLength',[0.02,0.02],'FontSize',38);
%axis off;
set(gca, 'FontName', 'Times');
hold on;
%%
figure;

[HH,LH,hh,lh] = comp_luttWF(Q.VB.WF);

p1=patch(isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(LH(:,:,:,3)).^2,0.05));
isonormals(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(LH(:,:,:,3)).^2, p1); 
daspect([1 1 1]);
set(p1, 'FaceColor', [0 0 1], 'EdgeColor', 'none');
axis tight;
%axis ([-10 10 -10 10 -60 60]);
camlight; lighting gouraud;
alpha(p1,0.5);
box on;
set(gca,'Linewidth',4,'TickLength',[0.02,0.02],'FontSize',38);
axis on;
hold on;
p2=patch(isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(HH(:,:,:,2)).^2,0.00005));
isonormals(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(HH(:,:,:,2)).^2, p2); 
daspect([2 2 1]);
set(p2, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
axis tight;
axis ([-20 20 -20 20 -20 20]);
camlight; lighting gouraud;
alpha(p2,0.5);
box on;
axis on;
hold on
set(gca, 'FontName', 'Times');
set(gca,'Linewidth',4,'TickLength',[0.02,0.02],'FontSize',38);