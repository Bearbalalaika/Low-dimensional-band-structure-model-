function varargout = ple(varargin)
% PLE M-file for ple.fig
% MB 22/06/2006
% ver 1.0
% call: [selection,title,label]=ple(listIn)
%      PLE, by itself, creates a new PLE or raises the existing
%      singleton*.
%      H = PLE returns the handle to a new PLE or the handle to
%      the existing singleton*.
%
%      PLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLE.M with the given input arguments.
%
%      PLE('Property','Value',...) creates a new PLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ple_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ple_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ple

% Last Modified by GUIDE v2.5 20-Jun-2006 12:51:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ple_OpeningFcn, ...
                   'gui_OutputFcn',  @ple_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ple is made visible.
function ple_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ple (see VARARGIN)

% Choose default command line output for ple
handles.output = -1;
handles.cuts=[];
handles.sel=[];
handles.clickedN={};
handles.list={};
if size(varargin)
    for k=1:size(varargin{1},2)
        handles.list{k}=char(varargin{1}(k));
    end;
else
    handles.list=evalin('base','dataFileNames');
end;
% Update handles structure
guidata(hObject, handles);
for k=1:size(handles.list,2)
    handles.accTime(1,k)=0;
end
% Update handles structure
guidata(hObject, handles);
set(handles.listbox1,'String',handles.list);
if evalin('base','exist(''dataLog'')');
    set(handles.uiLog,'String',evalin('base','dataLog'));
end
% UIWAIT makes ple wait for user response (see UIRESUME)
 uiwait(handles.figure1);

 
% --- Outputs from this function are returned to the command line.
function varargout = ple_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
varargout{1} = handles.output;
varargout{2} = get(handles.uiTitle,'String');
varargout{3} = get(handles.uiLabel,'String');
%varargout{4} = get(handles.uiNormalise,'Value');

if get(handles.uiNormalise,'Value')
    varargout{4} = [1 str2double(get(handles.uiNorm1,'String')) str2double(get(handles.uiNorm2,'String'))];
else
    varargout{4} = [0 str2double(get(handles.uiStep,'String')) 0];
end
delete(handles.figure1);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
handles.sel=get(hObject,'Value');
guidata(hObject, handles);
updateplot(hObject, eventdata, handles);
str=sprintf('%3.2f',handles.accTime(handles.sel(1,1)));
set(handles.uiAccTime,'String',str);



function str=guessVarName(str, n)
str=fliplr(char(str));
str=fliplr(str(1,strfind(str,'.')+1:end));
%if first char is a letter its ok - matlab did not modified it
%if it is a number, matlab added a X to a variable name
if length(str2num(str(1,1)))
    str=sprintf('X%s(:,%d)',str,n);%-------------------------    
else
    str=sprintf('%s(:,%d)',str,n);%-------------------------
end




function [xcut ycut]=getPLEcuts(strCut, handles)
%get cut energies from "uiCutEnergy" line by line - two values or one, and put into "cuts"
%TODO: use handels.cuts
%define ncuts based on handles 
xcut=[]; ycut=[];
ncuts=size(strCut,1);
%return if there is no cutting energy defined
if isequal(ncuts,0)
    return;
end;
for i=1:ncuts
    [str, str2]=strtok(strCut(i,:));
    cuts(i,1)=str2double(str);
    str2=strtok(str2);
    if size(str2,2)~=0
        cuts(i,2)=str2double(str2);
    else
        cuts(i,2)=0;
    end
end


%get current selection
sel=get(handles.listbox1,'Value');
if isequal(get(handles.uiReverse,'Value'),1)
    sel=fliplr(sel);
end

for i=1:size(sel,2)
    %get data
    n=sel(1,i);
    str=guessVarName(char(handles.list(1,n)),1);
    x=evalin('base',str);
    str=guessVarName(char(handles.list(1,n)),2);
    y=evalin('base',str);
    %switch X scale to eV or nm?
    if isequal(get(handles.uiEvNm,'Value'),1)
        x=1239.8./x; 
    end
    
    %get Y value at cut points
    %LAST WHEN NOT FOUND==SIZE
    for k=1:ncuts
        y(k,i)=0;
        %get index(es) of cutting points
        %x must be sorted !?
        in1=find(sort(x)<=cuts(k,1),1,'last');
        if isempty(in1)
            in1=0;
        end
        if cuts(k,2)==0
            %1 value -> get intensity at index
            if in1~=0 && in1~=size(x,1);
                yreg=y(in1,1);
                ycut(k,i)=sum(yreg);
            end
        else
            %2 values -> get region and calculate (integrated) intensity
            in2=find(sort(x)<=cuts(k,2),1,'last');
            if isempty(in2)
                in2=0;    
            end
            if in1==0 && in2==0
                yreg=0;
            elseif in1~=0 && in2~=0
                yreg=y(in1:in2,1);
            elseif in1==0
                yreg=y(1:in2,1);
            elseif in2==0
                yreg=y(in1:size(y,1),1);
            end 
            ycut(k,i)=sum(yreg);
        end   
    end %for k=1:ncuts
end

if ~exist('ycut')
    return;
end

xcut(1,:)=sel(1,:);
%if ple data was loaded transform sel to laser energy;
if evalin('base','exist(''dataLog'')') ;
    log= evalin('base','dataLog');
    total=size(handles.list,2)-size(log,2);
    %if 2 rows are the same width - probably a ple log
    if isequal(size(char(log(1,1)),2) , size(char(log(1,2)),2) )
        for i=1:size(ycut,2)
            if isequal(get(handles.uiReverse,'Value'),1)
                n=size(handles.list,2)-sel(1,i)+1;
            else
                n=sel(1,i);
            end
            xcut(1,i)=str2double(strtok(char(log(1,n))));
        end
        if isequal(get(handles.uiEvNm,'Value'),1)
            xcut=1239.8./xcut;
        end
    end

else
    xcut(1,:)=sel(1,:);
end



%--------------------------------------------------------
%--------------------------------------------------------
function [x y]=normalise(x,y,handles,i)
step=str2double(get(handles.uiStep,'String'));
stepN=str2double(get(handles.uiNormStep,'String'));
x1=str2double(get(handles.emin,'String'));
x2=str2double(get(handles.emax,'String'));
n1=str2double(get(handles.uiNorm1,'String'));
n2=str2double(get(handles.uiNorm2,'String'));
if isequal(n1,0)
    n1=x1;
end;
if isequal(n2,0)
    n2=x2;
end;
% change to eV or nm?
if isequal(get(handles.uiEvNm,'Value'),1)
   x=1239.8./x; 
end
%normalise to accumulation time?
if isequal(get(handles.uiUseAccTime,'Value'),1)
    if isequal(handles.accTime(1,n),0)
        %do nothing
    else
        %normalise to cps
        y=y/handles.accTime(1,n);
    end
end
%normalise to 1
if isequal(get(handles.uiNormalise,'Value'),1)
    %get index of x at n1,n2 energy
    in1=find(x>n1,1,'first');
    in2=find(x<n2,1,'last');   
    %get max and min values in the NORME, XRANGE of interest
    maxy=max(y(in1:in2,1));
    miny=min(y(intersect( find(x>x1),find(x<x2) ),1 ) );
    if i==0
        y=(y-miny)./(maxy-miny);
    else
        y=(y-miny)./(maxy-miny)+stepN*(i-1);
    end
else
    y=y+step*(i-1);
end;   



%--------------------------------------------------------
%--------------------------------------------------------
% --- Updates PL data plot
function updateplot(hObject, eventdata, handles)

%update PL plot (axes1)
set(gcf,'CurrentAxes',handles.axes1);
sel=get(handles.listbox1,'Value');
if isempty(sel)
    return;
end
%return if no selection
if isequal(sel,-1)
    return;
end

if isequal(get(handles.uiReverse,'Value'),1)
    sel=fliplr(sel);
end
%plot 1st selection as new plot
n=sel(1,1);
strx=guessVarName(char(handles.list(1,n)),1);
stry=guessVarName(char(handles.list(1,n)),2);
x=evalin('base',strx);
y=evalin('base',stry);
[x y]=normalise(x,y,handles,0);   
plot(handles.axes1,x,y );

if size(sel,2)>1
    %then hold and plot the rest
    hold on;
    for i=2:size(sel,2)
        n=sel(1,i);
        strx=guessVarName(char(handles.list(1,n)),1);
        stry=guessVarName(char(handles.list(1,n)),2);
        
        x=evalin('base',strx);
        y=evalin('base',stry);
        
        [x y]=normalise(x,y,handles,i);  
        plot(handles.axes1,x,y );
    end
    hold off;
end
%set limit values on axis from textboxes
lim=[str2double(get(handles.emin,'String')) str2double(get(handles.emax,'String'))];
set(handles.axes1,'XLim',lim);
lim=[str2double(get(handles.imin,'String')) str2double(get(handles.imax,'String'))];
set(handles.axes1,'YLim',lim);

%add lines indicating PLE cuts
% strCut=get(handles.uiCutEnergy,'String');
% ncuts=size(strCut,1);
% %return if there is no cutting energy defined
% if isequal(ncuts,0)
%     return;
% end;
% for i=1:ncuts
%     [str, str2]=strtok(strCut(i,:));
%     cuts(1,i)=str2double(str);
%     str2=strtok(str2);
%     if size(str2,2)~=0
%         cuts(2,i)=str2double(str2);
%     else
%         cuts(2,i)=0;
%     end
% end
% 
% disp(cuts);
% disp(lines);
% line(linesx,linesy)


%and update PLE plot...
updateplotPLE(hObject, eventdata, handles);



%--------------------------------------------------------
%--------------------------------------------------------
% --- Updates PLE data plot
function updateplotPLE(hObject, eventdata, handles)
%set as current: PLE plot (axes2)
%set(gcf,'CurrentAxes',handles.axes1);

strCut=get(handles.uiCutEnergy,'String');
ncuts=size(strCut,1);

[xcut ycut]=getPLEcuts(strCut, handles);

if ncuts~=0
    plot(handles.axes2,xcut(1,:),ycut,'o-');
end



% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor','white');



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);% The GUI is still in UIWAIT, us UIRESUME
else
    delete(handles.figure1);% The GUI is no longer waiting, just close it
end


% --- Executes on button press in uiAutoY.
function uiAutoY_Callback(hObject, eventdata, handles)
% hObject    handle to uiAutoY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axis 'auto y';
lim=get(handles.axes1,'YLim');
set(handles.imin,'String',num2str(lim(1,1),'%.6g')); 
set(handles.imax,'String',num2str(lim(1,2),'%.6g'));



function emin_Callback(hObject, eventdata, handles)
lim=[str2double(get(handles.emin,'String')) str2double(get(handles.emax,'String'))];
set(handles.axes1,'XLim',lim);

function emin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function emax_Callback(hObject, eventdata, handles)
lim=[str2double(get(handles.emin,'String')) str2double(get(handles.emax,'String'))];
set(handles.axes1,'XLim',lim);

function emax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function imax_Callback(hObject, eventdata, handles)
lim=[str2double(get(handles.imin,'String')) str2double(get(handles.imax,'String'))];
set(handles.axes1,'YLim',lim);

function imax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function imin_Callback(hObject, eventdata, handles)
lim=[str2double(get(handles.imin,'String')) str2double(get(handles.imax,'String'))];
set(handles.axes1,'YLim',lim);

function imin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btAutoScale.
function btAutoScale_Callback(hObject, eventdata, handles)
axis tight;
lim=get(handles.axes1,'XLim');
set(handles.emin,'String',num2str(lim(1,1),'%.6g') );
set(handles.emax,'String',num2str(lim(1,2),'%.6g') );
lim=get(handles.axes1,'YLim');
set(handles.imin,'String',num2str(lim(1,1),'%.6g') );
set(handles.imax,'String',num2str(lim(1,2),'%.6g') );


function uiTitle_Callback(hObject, eventdata, handles)
% hObject    handle to uiTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uiTitle as text
%        str2double(get(hObject,'String')) returns contents of uiTitle as a double


% --- Executes during object creation, after setting all properties.
function uiTitle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function uiLabel_Callback(hObject, eventdata, handles)
% hObject    handle to uiLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uiLabel as text
%        str2double(get(hObject,'String')) returns contents of uiLabel as a double


% --- Executes during object creation, after setting all properties.
function uiLabel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in uiNormalise.
function uiNormalise_Callback(hObject, eventdata, handles)
% hObject    handle to uiNormalise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of uiNormalise
updateplot(hObject, eventdata, handles);
axis 'auto y';
lim=get(handles.axes1,'YLim');
set(handles.imin,'String',num2str(lim(1,1),'%.4g') );
set(handles.imax,'String',num2str(lim(1,2),'%.4g') );

function uiNorm1_Callback(hObject, eventdata, handles)
% hObject    handle to uiNorm1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uiNorm1 as text
%        str2double(get(hObject,'String')) returns contents of uiNorm1 as a double
updateplot(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function uiNorm1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function uiNorm2_Callback(hObject, eventdata, handles)
% hObject    handle to uiNorm2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uiNorm2 as text
%        str2double(get(hObject,'String')) returns contents of uiNorm2 as a double
updateplot(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function uiNorm2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in uiOK.
function uiOK_Callback(hObject, eventdata, handles)
if isempty(handles.sel)
    handles.output=-1; %return -1 when no selection has been made
else
    handles.output=handles.sel;
end
guidata(hObject, handles);
uiresume(handles.figure1);



function uiStep_Callback(hObject, eventdata, handles)
% hObject    handle to uiStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uiStep as text
%        str2double(get(hObject,'String')) returns contents of uiStep as a double
updateplot(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function uiStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uiStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in uiMake.
function uiMake_Callback(hObject, eventdata, handles)

TITLE=char(get(handles.uiTitle,'String'));
LEFTTOP=char(get(handles.uiLabel,'String'));
NORMALISE=char(get(handles.uiNormalise,'Value'));     
STEPY=str2double(get(handles.uiStep,'String')); 
STEPYN=str2double(get(handles.uiNormStep,'String')); 
n1=str2double(get(handles.uiNorm1,'String'));   
n2=str2double(get(handles.uiNorm2,'String'));         
x1=str2double(get(handles.emin,'String')); 
x2=str2double(get(handles.emax,'String')); 
sel=handles.sel;
dataFileNames=handles.list;

%return if no selection
if isequal(sel,-1)
    return;
end
% 
if isequal(get(handles.uiReverse,'Value'),1)
    sel=fliplr(sel);
end
%make array of selected data names (get rid of filename extension)
for i=1:size(sel,2)
    n=sel(1,i);  
    str=fliplr(char(dataFileNames(1,n)));
    str=fliplr(str(1,strfind(str,'.')+1:end));
    selNames(1,i)=cellstr(str);
end;
%rearange data?
selNames=rearrange(selNames);
if isempty(selNames)
    return;
end;

%plot
fig = figure('Name',TITLE);
axes1 = axes('Parent',fig);
hold on;
for i=1:size(sel,2) 
    str=char(selNames(1,i));
    if length(str2num(str(1,1)))
            strx=sprintf('X%s(:,1)',str);%-------------------------
            stry=sprintf('X%s(:,2)',str);%-------------------------
        else
            strx=sprintf('%s(:,1)',str);%-------------------------
            stry=sprintf('%s(:,2)',str);%-------------------------
    end
    x=evalin('base',strx);
    y=evalin('base',stry);
    if isequal(get(handles.uiEvNm,'Value'),1)
        x=1239.8./x; 
    end
    %add time norm  
    if isequal(NORMALISE,1)
        %get index of CDX at NORME energy
        in1=find(x>n1,1,'first');
        in2=find(x<n2,1,'last');
        %get max and min values in the n1 n2, and x1 x2 range of interest
        maxy=max(y(in1:in2) );
        miny=min(y(intersect( find(x>x1),find(x<x2) ) ));
        %normalise y
        y=(y-miny)/(maxy-miny)+STEPYN*(i-1);
    else
        y=y+STEPY*(i-1);
    end;
    h=plot(axes1,x,y);  
end
hold off;

box on;
axis ([str2double(get(handles.emin,'String')) str2double(get(handles.emax,'String')) str2double(get(handles.imin,'String')) str2double(get(handles.imax,'String'))]);
% axis 'auto y';
xlabel('Energy (eV)');
if isequal(NORMALISE,1)
    ylabel('Normalised Intensity (arb.u.)');
else
    ylabel('Intensity (arb.u.)');
end
title(axes1,TITLE);
%legend(selNames);
annotation(fig,'textbox', [0.15 0.8 0.1 0.1], 'String',LEFTTOP, 'Margin', 0.001, 'FitHeightToText','on');
str=date;
annotation(fig,'textbox', [0.8 0.002 0.16 0.03], 'String',str, 'Margin', 0.001, 'FitHeightToText','on', 'LineStyle', 'none');
str=evalin('base','dataPathName');
str=strcat(str,' (', selNames(1,1),' - ', selNames(1,size(sel,2)),') ');
%annotation(fig,'textbox', [0.1 0.002 0.16 0.03], 'String',str, 'Margin', 0.001, 'FitHeightToText','on', 'LineStyle', 'none','Interpreter','none');



function uiAccTime_Callback(hObject, eventdata, handles)
% hObject    handle to uiAccTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uiAccTime as text
%        str2double(get(hObject,'String')) returns contents of uiAccTime as a double

%get sel from listbox
tmp=get(handles.listbox1,'Value');
sel=tmp(1,1);
handles.accTime(1,sel)=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
updateplot(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function uiAccTime_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in uiUseAccTime.
function uiUseAccTime_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of uiUseAccTime
updateplot(hObject, eventdata, handles);


function uiNormStep_Callback(hObject, eventdata, handles)
% hObject    handle to uiNormStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uiNormStep as text
%        str2double(get(hObject,'String')) returns contents of uiNormStep as a double
updateplot(hObject, eventdata, handles);

n=size(handles.sel,2);
step=str2double(get(hObject,'String'));
ymin=0;
ymax=ymin+step*(n+1);
ylim([ymin ymax]);

lim=get(handles.axes1,'YLim');
set(handles.imin,'String',num2str(lim(1,1),'%.4g') );
set(handles.imax,'String',num2str(lim(1,2),'%.4g') );


% --- Executes during object creation, after setting all properties.
function uiNormStep_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in uiReverse.
function uiReverse_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of uiReverse
sel=handles.sel;
list=handles.list;
list=fliplr(list);
handles.list=list;
sel=size(list,2)-sel+1;
handles.sel=sel;
% Update handles structure
guidata(hObject, handles);
set(handles.listbox1,'String',handles.list);
set(handles.listbox1,'Value',handles.sel);
updateplot(hObject, eventdata, handles);



% --- Executes on button press in uiEvNm.
function uiEvNm_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of uiEvNm
lim=[str2double(get(handles.emin,'String')) str2double(get(handles.emax,'String'))];
lim=1239.8./lim;
lim=fliplr(lim);
set(handles.axes1,'XLim',lim);
str=sprintf('%.5f',lim(1,1)); set(handles.emin,'String',str);
str=sprintf('%.5f',lim(1,2)); set(handles.emax,'String',str);


str=get(handles.uiCutEnergy,'String');

ncuts=size(str,1);
strOut='';txt='';
for i=1:ncuts
    [str1, str2]=strtok(str(i,:));
    cut1=str2double(str1);
    str2=strtok(str2);
    if size(str2,2)~=0
        cut2=str2double(strtok(str2));
        strOut=sprintf('%.4f %.4f',1239.8./cut1,1239.8./cut2);
    else
        strOut=sprintf('%.4f',1239.8./cut1 );
    end
    txt(i,:)=strOut;
end
set(handles.uiCutEnergy,'String',txt);
updateplot(hObject, eventdata, handles);



function uiCutEnergy_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of uiCutEnergy as text
%        str2double(get(hObject,'String')) returns contents of uiCutEnergy as a double

strCut=get(hObject,'String');
ncuts=size(strCut,1);
%return if there is no cutting energy defined

for i=1:ncuts
    [str, str2]=strtok(strCut(i,:));
    cuts(i,1)=str2double(str);
    str2=strtok(str2);
    if size(str2,2)~=0
        cuts(i,2)=str2double(str2);
    else
        cuts(i,2)=0;
    end
end
handles.cuts=cuts;
% Update handles structure
guidata(hObject, handles);


updateplot(hObject, eventdata, handles);
%updateplotPLE(hObject, eventdata, handles);




% --- Executes during object creation, after setting all properties.
function uiCutEnergy_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%--------------------------------------------------------------------------
%--------------------------  PLE SECTION  ---------------------------------
%--------------------------------------------------------------------------








% --- Executes on button press in uiMakePLE.
function uiMakePLE_Callback(hObject, eventdata, handles)
TITLE=char(get(handles.uiTitle,'String'));
LEFTTOP=char(get(handles.uiLabel,'String'));
NORMALISE=char(get(handles.uiNormalise,'Value'));     
STEPY=str2double(get(handles.uiStep,'String')); 
STEPYN=str2double(get(handles.uiNormStep,'String')); 
n1=str2double(get(handles.uiNorm1,'String'));   
n2=str2double(get(handles.uiNorm2,'String'));         
x1=str2double(get(handles.emin,'String')); 
x2=str2double(get(handles.emax,'String')); 
sel=handles.sel;
dataFileNames=handles.list;

strCut=get(handles.uiCutEnergy,'String');
ncuts=size(strCut,1);
[xcut ycut]=getPLEcuts(strCut, handles);

%plot
fig = figure('Name',TITLE);
axes1 = axes('Parent',fig);
hold on;
    h=plot(axes1,xcut(1,:),ycut,'o-');  
hold off;
box on;
axis 'auto';
xlabel('Excitation Energy (eV)');
ylabel('PLE Intensity(arb.u.)');
title(axes1,TITLE);
legend(strCut);
annotation(fig,'textbox', [0.15 0.8 0.1 0.1], 'String',LEFTTOP, 'Margin', 0.001, 'FitHeightToText','on');
str=date;
annotation(fig,'textbox', [0.8 0.002 0.16 0.03], 'String',str, 'Margin', 0.001, 'FitHeightToText','on', 'LineStyle', 'none');
str=evalin('base','dataPathName');
str=strcat(str,' (', selNames(1,1),' - ', selNames(1,size(sel,2)),') ');
annotation(fig,'textbox', [0.1 0.002 0.16 0.03], 'String',str, 'Margin', 0.001, 'FitHeightToText','on', 'LineStyle', 'none','Interpreter','none');



%--------------------------------------------------------------------------
%------------------------  RECENTLY ADDED  --------------------------------
%--------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function uiLogPath_CreateFcn(hObject, eventdata, handles)
str=evalin('base','dataPathName');
set(hObject,'String',str);


% --- Executes on button press in uiExportPLE.
function uiExportPLE_Callback(hObject, eventdata, handles)

sel=handles.sel;
dataFileNames=handles.list;

strCut=get(handles.uiCutEnergy,'String');
ncuts=size(strCut,1);
[xcut ycut]=getPLEcuts(strCut, handles);

OUT=[xcut; ycut]';
str=inputdlg('Variable name:','Input',1,{'plexy'});
if isempty(str)
    return
end
assignin('base',char(str),OUT);
disp('ok');

