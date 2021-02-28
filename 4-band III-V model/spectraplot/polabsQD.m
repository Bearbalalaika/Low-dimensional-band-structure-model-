%clear all;
%load UNEQUAL_QD1_9nm_QD2_0d5nm_DQD_4dz_Barrier_4nm.mat;
Q = exciton_pert(Q);

[mlh, Elh] = absQD_II(Q.CB.WF, Q.VB.WF, [0 0 1], Q.CB.E, Q.VB.E);

%
Eb = 0;
if isfield(Q, 'X')
    disp('Exciton binding included')
    Eb = Q.X.Ex; 
end
    
Elh = Elh-Eb;
%}

[mhh, Ehh] = absQD_II(Q.CB.WF, Q.VB.WF, [1 1 0], Q.CB.E, Q.VB.E);

%
Eb = 0;
if isfield(Q, 'X')
    disp('Exciton binding included')
    Eb = Q.X.Ex; 
end
    
Ehh = Ehh-Eb;

%}

Elh=Elh./1000;
Ehh=Ehh./1000;
maxlh=max(max(mlh));




figure; stem(1.5192+Elh, mlh/maxlh, 'bo','linewidth',2);
hold on; stem(1.5192+Ehh, mhh/maxlh, 'ro','linewidth',2);


set(gca,'fontsize',20);

set(get(gca,'YLabel'),'String','I [a.u.]','Color','k','fontsize',24);
set(get(gca,'XLabel'),'String','Ev [eV]','fontsize',24);

hold on;

for i=1:5
    for j=1:2:24
      %  stem(1.5192+Elh(i,j), mlh(i,j), 'b');
      t=(j+1)/2;
        if (mhh(i, j)>0.008) h = text(1.519+Elh(i, j), n+1.1*mhh(i, j), sprintf('e_{%d}h_{%d}', i, t)); set(h, 'color', [1 0 0]*0.75, 'fontsize', 12, 'fontweight', 'bold');end
        if (mlh(i, j)>0.01) h = text(1.519+Ehh(i, j), n+1.1*mlh(i, j), sprintf('e_{%d}h_{%d}', i, t)); set(h, 'color', [0 0 1]*0.75, 'fontsize', 12, 'fontweight', 'bold');end

    end
end

%light and heavy hole spectra part

