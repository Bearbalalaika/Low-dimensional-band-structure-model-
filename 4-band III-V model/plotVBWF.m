%clear all;
%load UNEQUAL_7nm_DQD_20dz_Barrier.mat;
[HH,LH,hh,lh] = comp_luttWF(Q.VB.WF);
figure;


p1=patch(isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(LH(:,:,:,1)).^2, 0.5));
isonormals(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(LH(:,:,:,1)).^2, p1); 
daspect([1 1 1]);
set(p1, 'FaceColor', [0 0 1], 'EdgeColor', 'none');
%axis tight;
%axis ([-10 10 -10 10 -80 80]);
camlight; lighting gouraud;
alpha(p1,1);
box on;
%max(max(max(abs(HH(:,:,:,1)))).^2)/2.72
axis on;
hold on;
p2=patch(isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(HH(:,:,:,1)).^2, max(max(max(abs(HH(:,:,:,1)))).^2)/2.72));
isonormals(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(HH(:,:,:,1)).^2, p2); 
daspect([1 1 1]);
set(p2, 'FaceColor', [1 0 0], 'EdgeColor', 'none');
axis tight;
axis ([-10 10 -10 10 -40 40]);
camlight; lighting gouraud;
alpha(p2,1);
box on;
axis on;



set(gca,'Fontsize',20,'Linewidth',2);
set(gca, 'FontName', 'Times');
set(gca,'FontSize',40);

%%
figure;
p0=patch(isosurface(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(Q.CB.WF(:,:,:,1)).^2, max(max(max(abs(Q.CB.WF(:,:,:,1)))).^2)/2.72));
isonormals(Q.lrpot.x, Q.lrpot.y, Q.lrpot.z, abs(Q.CB.WF(:,:,:,1)).^2, p0); 
daspect([1 1 1]);
set(p0, 'FaceColor', [0 1 0], 'EdgeColor', 'none');
%axis tight;
%axis ([-15 15 -15 15 -80 80]);
camlight; 
lighting gouraud;
%lighting phong;
%lighting flat;
alpha(p0,1);box on;
%axis on;
axis ([-10 10 -10 10 -230 230]);
%set(gca, 'XTick', [],'YTick', []);
%set(gca, 'ZTickLabel','','Linewidth',1.5,'TickLength',[0.02,0.02]);
set(gca,'Linewidth',1.5,'TickLength',[0.02,0.02],'FontSize',16);
set(gca, 'FontName', 'Times');
set(gca,'FontSize',40);
%axis off;