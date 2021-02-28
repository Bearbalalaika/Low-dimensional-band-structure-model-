%clear all;
%load UNEQUAL_QD1_9nm_QD2_0d5nm_DQD_4dz_Barrier_4nm.mat;
%Q = exciton_pert(Q);

[Elh, Ilh, Enlh, Mnlh] = absII_cont(Q, [0 0 1]);

[Ehh, Ihh, Enhh, Mnhh] = absII_cont(Q, [1 1 0]);


Elh=Elh./1000;
Ehh=Ehh./1000;

%1.5192
figure; stem(1.519+Elh, Ilh, 'bo','linewidth',2);
hold on; stem(1.519+Ehh, Ihh, 'ro','linewidth',2);

set(gca,'fontsize',24);
xlim([1.585 1.62]);

set(get(gca,'YLabel'),'String','I [a.u.]','Color','k','fontsize',28);
set(get(gca,'XLabel'),'String','Ev','fontsize',28);


%light and heavy hole spectra part

