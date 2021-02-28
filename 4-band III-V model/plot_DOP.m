function  plot_DOP( specV,specH ,varargin)
% plot_DOP  This function plots a spectrum linearly resolved in polarization 
% and it's degree of polarization.
% OPTIONS:  -'PlotLim' , type->array: give plot limits as function arguements
%           -'Font' , type->char: specify text font as function argument
%           -'FontSize', type->double: give text fontsize as function argument
%           -'AxisPos', type->array: give axis position as function
%           argument
%           -'DisplayQ', type->array(x0,Q): display quality factor Q at position
%           x0
%           -'LogInfo', type->array(temp,power): give log info numbers for 
%           temperature and power as function arguments 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update 10/02/2014 (Clement): added doc for function arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=299792458;
hplank=4.13566751691*1e-15;
% Choose energy lower and upper limit
minE=min(specV.E);
maxE=max(specV.E);

% Creation of default values
fontname='Helvetica';
F1=30;
axispos= [100,500,1200,700];
Lim=0;
Q=0;
num_temp=10;
num_pow=8;
% Reads arguemnts list and specifies new values
for i=1:length(varargin)
   if strcmp(varargin(i),'PlotLim')
      Lim=varargin{i+1}; 
      minE=Lim(1); maxE=Lim(2);
   elseif  strcmp(varargin(i),'Font') 
      fontname=varargin{i+1};
   elseif  strcmp(varargin(i),'FontSize') 
      F1=varargin{i+1};
   elseif  strcmp(varargin(i),'AxisPos') 
      axispos=varargin{i+1};
   elseif  strcmp(varargin(i),'DisplayQ') 
      temp=varargin{i+1}; 
      x0=temp(1); Q=temp(2);
   elseif  strcmp(varargin(i),'LogInfo') 
      temp=varargin{i+1}; 
      num_temp=temp(1); num_pow=temp(2);
   end
end

% removal of noise for plot
IV=(specV.I);
IH=(specH.I);
noise=mean(IH(end-10:end)); %can be modified 
IV=IV-noise;
IH=IH-noise;
E=specV.E;
nmax=max(max(IV),max(IH));
nmin=min(min(IV),min(IH));

set(0,'defaultaxesfontname',fontname);

Line=3;  % default frame width
lw=3.2;  %default linewidth
h=figure;
set(h, 'position', axispos);


% plot of V polarized spectrum and dgree of polarization (DOP)
[AX,H1,H2]=plotyy(E,IV,E,(IV-IH)./(IV+IH+20),'plot');
% axis handle for V and H polarization
ax1=AX(1);
% axis handle for DOP
ax2=AX(2);
% plot of H polarized spectrum
H3=line(E,IH,'Color','k','linewidth',lw,'parent',ax1);

uistack(H1,'top');
uistack(H3,'top');

% Properties of axis for V and H spectra
set(ax1,'XColor','k','YColor','k','Box','off')
set(ax1,'TickLength',[0.022 0.0045],'YMinorTick','on','TickDir','in')
set(ax1,'Xlim',[minE maxE]);
%axis energy to wavelength

lab=(get(ax1,'XTick'));
set(ax1, 'XTickLabel',num2str((c*hplank)./lab'*1e9,3));
set(ax1,'XTickMode','manual');

tick=linspace(0,nmax,4);
set(ax1,'YTick',tick,'YTickLabel',[]);
set(ax1,'Ylim',[-0.1*nmax, nmax+0.5*nmax]);
set(get(ax1,'Ylabel'),'String','I [a.u.]','fontsize',F1,'FontName',fontname)
set(get(ax1,'Xlabel'),'String','Wavelength [nm]','fontsize',F1,'FontName',fontname)
set(ax1,'FontSize',F1,'XAxisLocation','Top','Linewidth',Line);
% Porperties of V polarized graphic object
set(H1,'color','r');
set(H1,'linewidth',lw);

% Properties of axis for DOP
set(ax2,'YColor',[0.4 0.4 0.4]);
set(ax2,'fontsize',F1);
set(ax2,'Xlim',[minE maxE]);
set(ax2,'YTick',[-1 0 1]);
set(ax2,'Ylim',[-8  1]);
set(get(ax2,'YLabel'),'String','DOP','fontsize',F1,'Position',get(get(ax2,'YLabel'),'Position') + [-0.000 3.5 -0.5],'FontName',fontname)
set(get(ax2,'Xlabel'),'String','Energy [eV]','fontsize',F1,'FontName',fontname)
set(ax2,'TickLength',[0.022 0.0045],'Linewidth',Line,'TickDir','in')

% Plot of line showing zero DOP
hl2 = line(specV.E,zeros(1,length(specV.E)),'Color',[0.5 0.5 0.5],'linestyle','--','linewidth',1.1,'parent',ax2);
hl3 = line(specV.E,-1*ones(1,length(specV.E)),'Color',[0.5 0.5 0.5],'linestyle','-','linewidth',1.1,'parent',ax2);

% Porperties of DOP (graphic object)
set(H2,'color',[0.3 0.3 0.3]);
set(H2,'linewidth',lw);
uistack(H2,'top');

% Legend and additional text on figure
leg=legend(ax1,'I_V','I_H');
legtxt=findobj(leg,'type','text');
set(leg, 'Position', [0.7728  ,  0.5105  ,  0.1241  ,  0.1372],'FontName',fontname);
set(legtxt(2),'color','r')
legend(ax1,'boxoff');
eval( [ 'text(maxE-0.2*(maxE-minE),nmax*0.5,specV.info',num2str(num_temp),',''FontSize'',F1,''FontName'',fontname);'] );
eval( [ 'text(maxE-0.2*(maxE-minE),nmax*0.35,specV.info',num2str(num_pow),',''FontSize'',F1,''FontName'',fontname);'] );
set(ax1, 'ActivePositionProperty' , 'OuterPosition','Units','Normalized');
set(ax2, 'ActivePositionProperty' , 'OuterPosition','Units','Normalized');

% Corrects scalling of axes
   scalefactor = .98;

g = get(ax1,'Position');
g(1:2) = g(1:2) + (1-scalefactor)/2*g(3:4);
g(3:4) = scalefactor*g(3:4);
set(ax1,'Position',g);

g = get(ax2,'Position');
g(1:2) = g(1:2) + (1-scalefactor)/2*g(3:4);
g(3:4) = scalefactor*g(3:4);
set(ax2,'Position',g);


set(gca,'Position',[0.07 0.12 0.83 0.76]);
%
(c*hplank)./lab'*1e9;

if(Q~=0)
    text(x0+0.01*(maxE-minE),nmax*0.9,['Q=',num2str(Q,4)],'FontSize',20,'FontName',fontname);
end

end

