
% This script plots matrix elements of dipol matrix which correspond to HH
% and LH tansitions
%Q = exciton_pert(Q);
[ml,mhh, E] = LH_HH_spectra(Q.CB.WF, Q.VB.WF,  Q.CB.E, Q.VB.E);

E=E./1000;

figure;
n=0;
%maxml=max(max(mhh));
maxml=max(max(mhh));
stem(1.5192+E, ml/maxml, 'bo','linewidth',2);
hold on; stem(1.5192+E, mhh/maxml, 'ro','linewidth',2);
set(gca,'fontsize',28);
hold on;
%h = getframe;
LtrNf=zeros(20);
HtrNf=zeros(20);
for i=1:size(E, 1)
    for j=1:2:size(E, 2)
        tr=(j+1)/2;
        
        Ltr(i, tr)=ml(i, j+1);
        Htr(i, tr)=mhh(i, j+1);
        %if (ml(i, j)>0.01) h = text(1.519+E(i, j), n+1.45*ml(i, j+1), sprintf('e_{%d}h_{%d}', i, 1+(j-1)*0.5)) ; set(h, 'color','blue','fontsize', 18, 'fontweight', 'bold'); end
        %if (mhh(i, j)>0.01) h = text(1.519+E(i, j), n+1.5*mhh(i, j+1), sprintf('e_{%d}h_{%d}', i, 1+(j-1)*0.5)) ; set(h, 'color','red', 'fontsize', 18, 'fontweight', 'bold'); end
        LtrN(i, tr)=Ltr(i, tr)/maxml;
        HtrN(i, tr)=Htr(i, tr)/maxml;
        
        if ((i<20) & (tr<20)) LtrNf(i, tr)=LtrN(i, tr);  HtrNf(i, tr)=HtrN(i, tr); end
    
    

   
        
    end
    
    
end
set(get(gca,'YLabel'),'String','Intensity [a.u.]','Color','k','fontsize',28);
set(get(gca,'XLabel'),'String','Energy [eV]','fontsize',28);
axis([1.63 1.65 -0.05 1.1]);
set(gca,'Linewidth',2,'FontSize',28);
%axis off;
set(gca, 'FontName', 'Times');





figure;
bar3(LtrN,'b')
set(gca,'Linewidth',2,'FontSize',24);
hold on;
%axis off;
set(gca, 'FontName', 'Times');
bar3(HtrN,'r')
%mhh
