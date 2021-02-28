function out = absII_cont_pol(Q, n)
if nargin<2
    n = 0;
end

    [E, Iz, En, Mz] = absII_cont(Q, [0 0 1]);
%    [E, Iy] = absII_cont(Q, [0 1 0]);
    [E, Ix, En, Mx] = absII_cont(Q, [1 0 0]);


%Nf = max([Ix Iy Iz]);
Nf = max([Ix Iz]);

if (nargin<2) figure; end

hold on;
        


%Ixy = (Ix+Iy)/2;
%P = (Ixy - Iz)./(Ixy + Iz);
%h4 = plot(1519+E, P, 'g'); 
%set(h4, 'linewidth', 4);
E=E./1000;%%modified 0506

%h1 = plot(1.519+E, n + Iz/Nf, 'b','LineWidth',2);
%h2 = plot(1.519+E, n + Ix/Nf, 'r','LineWidth',2);

%h1 = plot(1519+E, n + Iz/Nf, 'b');
%h2 = plot(1519+E, n + Ix/Nf, 'r');

out(1, :) = 1.519+E;%%modified 0506
out(2, :) = Ix/Nf;
out(3,  :) = Iz/Nf;




for i=1:length(out(3,  :))

check=out(3,  i) + out(2, i);

if check == 0
    out(4,  i) = 0;
    
else
    
    out(4,  i) =(out(3,i)- out(2,i))./(out(3,  i)+out(2, i)+0.03);
    
end

end

% plot of V and H polarized spectrum
subplot(2,1,1);
[AX,h1,h2]=plotyy(1.519+E, n + Iz/Nf,1.519+E,n + Ix/Nf,'plot')
%1.519

axis(AX,[1.59 1.62 -0.05 1.1]);
set(gca, 'FontName', 'Times');
AX=gca;
%xlabel(AX,'Ev');
%ylabel(AX,'I [a.u.]');

%subplot(2,1,1);
%h3 = plot(1.519+E,out(4,:),'k','LineWidth',2)


size(E)



set(get(AX,'YLabel'),'String','Intensity [a.u.]','Color','k','fontsize',28);
set(get(AX,'XLabel'),'String','Energy [eV]','fontsize',28);
set(gca,'linewidth',3,'fontsize',24);

%h3 = plot(1519+E, 0 + Iy/Nf, 'b:'); 
%subplot(2,1,2);
set(h1, 'linewidth', 2,'Color','b');
set(h2, 'linewidth', 2,'Color','r');



%axis([1580 1615 -1 3.25]); box on;
%axis([1580 1615 -0.1 1.5]); box on;

sx = size(En, 2);
En = En(:, 1:2:sx);
Mx = Mx(:, 1:2:sx);
Mz = Mz(:, 1:2:sx);



%{
for i=1:size(En, 1)
    for j=1:size(En, 2)
        if (Mx(i, j)>0.0000008) h = text(1.519+En(i, j)/1000, n+400.6*Mx(i, j), sprintf('e_{%d}h_{%d}', i, j)); set(h, 'color', [1 0 0]*0.75, 'fontsize', 12, 'fontweight', 'bold');end
        if (Mz(i, j)>0.0000001) h = text(1.519+En(i, j)/1000, n+300.7*Mz(i, j), sprintf('e_{%d}h_{%d}', i, j)); set(h, 'color', [0 0 1]*0.75, 'fontsize', 12, 'fontweight', 'bold');end
    end
end

%}

% plot of dgree of polarization (DOP)
subplot(2,1,2);


plot(1.519+E, out(4,:),'k','LineWidth',2)

axis([1.59 1.62 -1 1]);
xlabel('Energy [eV]','fontsize',18)
ylabel('DOLP','fontsize',18)

set(gca, 'FontName', 'Times');
%set(h3, 'linewidth', 2);

%axis(ax2,[1.550 1.660 -1 1]);
%ax2 = gca;


%set(get(ax2,'YLabel'),'String','DOP','fontsize',18)
%set(get(ax2,'XLabel'),'String','Ev','fontsize',18)


%set(get(ax2,'Xlabel'),'String','Energy [eV]','fontsize',10)




    