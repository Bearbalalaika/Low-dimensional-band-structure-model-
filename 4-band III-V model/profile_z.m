

figure;

profile = zeros(size(z));

for i=1:length(z)
    profile(i)=Q.lrpot.C(15, 15,i);
end
%figure;
plot (z,profile);
set(get(gca,'YLabel'),'String','Aluminium content','Color','k','fontsize',32);
set(get(gca,'XLabel'),'String','Z[nm]','fontsize',32);
set(gca,'fontsize',26);
hold on;
set(gca, 'FontName', 'Times');
%%
% CB potential

figure;

profile = zeros(size(z));
energyprof=zeros(size(z));
for i=1:length(z)
    profile(i)=Q.lrpot.C(25, 25,i);
   % energyprof(i)=(1.36*profile(i)+0.22*profile(i)^2)*1000 ;
   %energyprof(i)= (polyval([0.22 1.36 0], profile(i)))*1000;
   energyprof(i)= (polyval(Q.mparam.EgAlGaAs, profile(i))+polyval(Q.mparam.VBOAlGaAs,profile(i)))*1000;
    %energyprof(i)=(1.36*0.047+0.22*0.04^27)*1000 ;
   
    
end
%figure;

CBenergyprof=energyprof+ Efieldz*z*1000;

plot (z,CBenergyprof);
set(get(gca,'YLabel'),'String','Energy [meV]','Color','k','fontsize',32);
set(get(gca,'XLabel'),'String','Z[nm]','fontsize',32);
set(gca,'fontsize',26);
hold on;
set(gca, 'FontName', 'Times');

%%

% VB potential
figure;

profile = zeros(size(z));
energyprof=zeros(size(z));
for i=1:length(z)
    profile(i)=Q.lrpot.C(20, 20,i);
    %energyprof(i)=(1.36*profile(i)+0.22*profile(i)^2)*1000 ;
    %energyprof(i)=polyval([0.22 1.36 0],profile(i) )*1000;
    energyprof(i)=polyval(Q.mparam.VBOAlGaAs,profile(i) )*1000;
end
%figure;

VBenergyprof=energyprof + Efieldz*z*1000;

plot (z,VBenergyprof);
set(get(gca,'YLabel'),'String','Energy [meV]','Color','k','fontsize',32);
set(get(gca,'XLabel'),'String','Z[nm]','fontsize',32);
set(gca,'fontsize',26);
hold on;
set(gca, 'FontName', 'Times');


%%



%%
profile = zeros(length(x),length(z));

for i=1:length(x)
    

for j=1:length(z)
    profile(i,j)=Q.lrpot.C(i, 40,j);
end
end

figure;
surf (z,x,profile);

set(gca,'Linewidth',1.5,'FontSize',28);
colorbar;
set(gca,'Linewidth',1.5,'FontSize',28);
%}
hold on;
set(gca, 'FontName', 'Times');
%%
%XY

profile = zeros(length(x),length(y));

for i=1:length(x)
    

for j=1:length(y)
    profile(i,j)=Q.lrpot.C(i, j,25);
end
end

figure;
surf (x,y,profile);

set(gca,'Linewidth',1.5,'FontSize',28);
colorbar;
set(gca,'Linewidth',1.5,'FontSize',28);
%}
set(gca, 'FontName', 'Times');