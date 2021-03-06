function varargout = ple2(varargin)
%
% Edit the above text to modify the response to help ple2

% Last Modified by GUIDE v2.5 26-Sep-2007 12:54:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ple2_OpeningFcn, ...
    'gui_OutputFcn',  @ple2_OutputFcn, ...
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


%--------------------------------------------------------------------------
% --- Executes just before ple2 is made visible.
function ple2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ple2 (see VARARGIN)

movegui(hObject,'northwest');

handles.nFiles = 0;             %number of items
handles.file_names    = [];     %loaded filenames..
handles.selectedLine  = 0;      %selected line label (ie. X-) (+2 added in update_ bcs. header)
handles.selectedIndex = 0;      %selected file in FileBox (+2 added in update_ to mark correct on bcs. headers)
handles.filteredIndex = [];     %filtered index F- get real index R by: R=h.i_f(F)

%fileDataType = struct('name','','x',[],'y',[]); PRELOADING IS NOT MUCH FASTER
%handles.allData = fileDataType;
%handles.alldata = [];           %data from directory
handles.data = [];              %data of dsplayed files ( pts(2048) x 2 x num spec)
handles.eng  = [];
handles.pow  = [];
handles.tac  = [];              % accumulation time

handles.tool = 'none';
handles.hSpec= [];              % = handles to lines;
handles.yStep = 0;              % step between spectra

handles.hMarkers  = [];
handles.markerPos = [];
handles.xroi = [];
handles.cutEng = 0;             % updated when PLE plot is drawn update_plot2()


set(hObject,'CurrentAxes',handles.Axes1);
handles.hMarkers(1,1) = vline(0,'b');
handles.hMarkers(1,2) = hline(0,'b');
set(handles.hMarkers,'Visible','off','Tag','Marker');

handles.output = hObject;       % Choose default command line output for ple2
guidata(hObject, handles);      % Update handles structure

%get dir of point_reader, append hgtable and add to matlab path
path(path,fullfile(fileparts(which('point_reader')),'hgtable')) ;

%load all data

list_dir(hObject,pwd);
filter_apply(hObject);

% UIWAIT makes ple2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%--------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = ple2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%% - small helpers, constants
function c = COLORPLOT(str)
c = [0.15 0.15 0.15];



%--------------------------------------------------------------------------
function Mask_Callback(hObject, eventdata, handles)


%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function Mask_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%% - FileBox -

%--------------------------------------------------------------------------
% --- Executes on selection change in FileBox.
function FileBox_Callback(hObject, eventdata, handles)
hFig = ancestor(hObject,'figure','toplevel');
handles = guidata(hFig);

items          = handles.nFiles;
oldSelection   = handles.selectedIndex;

index_selected = get(handles.FileBox,'Value') - 2;  % header is 2 lines
if index_selected > items                           % check selection
    return
end
index_selected = index_selected( index_selected > 0); % get only positie selection , to avoid selecting header
if isempty(index_selected)
    index_selected = oldSelection;
    set(hObject,'Value', oldSelection+2);
else
    handles.selectedIndex = index_selected;
    guidata(hObject, handles);                          % Update handles structure
end

% load data
index = handles.filteredIndex;
file  = handles.selectedIndex;
if min(file)<=0
    return
end
ifile=index(file);
filename = handles.file_names(ifile);
N = length(ifile);

for i=1:N
%       was: load data selected in filebox, then plot them
%       now: data is loadad at start to handles.allData(name,x,y)
%            and selected files are loaded from memory
%       Reverted: preloading does not improve plotting time
    data=load(char(filename{i}));   % load data
    data(:,1) = data(:,1) * 1000;           % change to meV
    y = data(:,2);
    handles.data(:,:,i)=data;       % matrix XY,XY,...
    %handles.data(:,:,i) = [handles.allData(ifile(i)).x handles.allData(ifile(i)).y];
end

set(handles.SelectionBox,'String',char(filename)); %list selected names
set(handles.SelectionBox,'Value',1);
guidata(hObject, handles);                          % Update handles structure
update_plot(hObject);
update_plot2(hFig);


%--------------------------------------------------------------------------
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over FileBox.
function FileBox_ButtonDownFcn(hObject, eventdata, handles)
hFig = ancestor(hObject,'figure','toplevel');
handles = guidata(hFig);
data = handles.data;
index = handles.filteredIndex;
file  = handles.selectedIndex;
if min(file)<=0
    return
end
ifile=index(file);
filename = handles.file_names(ifile);
%rename line
if strcmp(get(handles.figure1,'SelectionType'),'open')
    def = {sprintf('%g',handles.pow(ifile)) , sprintf('%g',handles.eng(ifile)) };
    answer = inputdlg({'Enter Power (uW):','Enter Lexc (nm):'},char(filename),1,def);
    if ~isempty(char(answer{1}))
        handles.pow(ifile) = str2num(answer{1});
    end
    if ~isempty(char(answer{2}))
        handles.eng(ifile)=str2num(answer{2});
    end
end
% Update handles structure
guidata(hObject, handles);
update_FileBox(hObject);



%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function FileBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%--------------------------------------------------------------------------
function update_FileBox(hObject)
%get file list form current_dir and mask
handles = guidata(hObject);
file   = handles.selectedIndex;
index  = handles.filteredIndex;
names  = cell(size(index));
eng    = cell(size(index));
pow    = cell(size(index));
tac    = cell(size(index));

num_files = size(handles.file_names,1);
if num_files<1
    set(handles.FileBox,'String','No files')
    return
end

if file < 1
    file = 1;
end

for i=1:size(index,1)
    ci = index(i);
    names{i} = char(handles.file_names(ci,1));
    eng{i}   = sprintf('%g', handles.eng(ci,1));
    pow{i}   = sprintf('%g', handles.pow(ci,1));
    tac{i}   = sprintf('%g', handles.tac(ci,1));
end

%from tablesetup:
set(handles.FileBox,'Value',file+2); %2 lines of headers
set(handles.FileBox,'FontName','FixedWidth','min',0,'max',2); % multiselection
formattable(handles.FileBox,{'File','P(uW)','Lex(nm)', 't(s)'}, names,pow,eng,tac);
set(handles.FileBox,'Value',file+2); %2 lines of headers
% Update handles structure
guidata(hObject, handles);



%--------------------------------------------------------------------------
function list_dir(hObject,dir_path)
handles = guidata(hObject);
mask=get(handles.Mask,'String');
mask=strcat('*',mask);

dir_struct = dir(fullfile(dir_path,mask));
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.file_names = sorted_names;
handles.filteredIndex = [1:size(handles.file_names,1)]';

%reset fields matrix for newly listed files
% PRELOADING waitbar
% h=waitbar(0,'Clearing...');  
% movegui(h,'northwest');
% set( get( get(h,'Children'),'Title') ,'Interpreter','none');

N=size(handles.file_names,1);
for i=1:N
    handles.pow(i,1) = 0;
    handles.eng(i,1) = 0;
    handles.tac(i,1) = 0;
%     following PRELOADING does not speed up drawing
%     filename = handles.file_names{i};
%     waitbar(i/N,h,char(filename));
%     data=load(char(filename));   % load data
%     data(:,1) = data(:,1) * 1000;           % change to meV
%     handles.allData(i).name = char(filename); 
%     handles.allData(i).x = data(:,1);
%     handles.allData(i).y = data(:,2);
end
%close(h); % PRELOADING waitbar

set(handles.FileBox,'Value',1);     % make default selection =1 in box (+2 in update_FileBox() )
handles.selectedIndex=-1;           % 2 lines of headers
guidata(hObject, handles);          % Update handles structure
update_FileBox(hObject);
contents = get(handles.FileBox,'String');%count number of items
items = size(contents,1);
handles.nFiles=items;
guidata(hObject, handles);          % Update handles structure


%--------------------------------------------------------------------------
function filter_apply(hObject)
handles = guidata(hObject);

names = handles.file_names;
index = handles.filteredIndex;
file  = handles.selectedIndex;
%filterpol = get(handles.Filter_Pol,'Value');
%filterpow = get(handles.Filter_Pow,'Value');
if min(file) <= 0
    file = 1;
end
if isempty(index) || isempty(names)
    return
end
num_files = size(names,1);

ifile = index(file);
name  = names(ifile);
%pow = handles.pow(ifile);
%eng = handles.eng(ifile);
index_eng=[];
index_pow=[];

index=1:num_files;

% if filterpol
%     index_pol = find(handles.polarization == pol);
%     index = intersect(index,index_pol);
% end
%
% if filterpow
%     index_pow = find(handles.power == pow);
%     index = intersect(index,index_pow);
% end

handles.filteredIndex= index';

%look for name and re-select it
for i=1:size(index,2)
    if strcmp(name,names(index(i)))
        break;
    end
end
handles.selectedIndex = i;
% Update handles structure
guidata(hObject, handles);

update_FileBox(hObject);








%% - Plot -


%--------------------------------------------------------------------------
function update_plot(hObject)
% Item selected in list box

hFig = ancestor(hObject,'figure','toplevel');
handles = guidata(hFig);
data = handles.data;
index = handles.filteredIndex;
file  = handles.selectedIndex;
yStep = handles.yStep;
if min(file)<=0
    return
end
ifile=index(file);
filename = handles.file_names(ifile);
N = length(ifile);

%set(handles.SelectionBox,'String',char(filename)); %list selected names


%set CURRENT AXES to axes1 - main plot
set(handles.figure1,'CurrentAxes',handles.Axes1);
set(handles.Axes1,'NextPlot','replacechildren');

rev = get(handles.PlotReversed,'Value');
nrm = get(handles.Normalize,   'Value');
bgd = get(handles.RemoveBgnd,  'Value');
tac = get(handles.ConsiderTac, 'Value');

for i=1:N
    if rev
        j = N - i + 1;
    else
        j = i;
    end
    x = data(:,1,j);            % change to meV
    y = data(:,2,j);
    if tac 
        y= y/handles.tac(ifile(j));
    end
    
    s = size(x,1);
    xi = x(1,1) : (x(s,1)-x(1,1))/(s.*2) : x(s,1);
    yi = interp1(x,y,xi);       % interpolate spetrum ( -> 2x points )

    n = xlim(handles.Axes1);
    in1=find(xi>n(1),1,'first');
    in2=find(xi<n(2),1,'last');
    if bgd
        miny = min(yi); % alll yi and not: ..(intersect( find(xi>n(1)),find(xi<n(2)) ) ));
        yi = yi - miny;
        handles.bgd = miny;
    end

    if nrm
        maxy = max( yi(in1:in2) );
        yi=(yi)/(maxy);
        yStep=handles.yStep/1000;
    end

    plot(xi,yi+yStep*(i-1),'Color',COLORPLOT);                %update plot in axes1
    if i==1
        set(handles.Axes1,'NextPlot','add');
    end
    hline=findobj(gca,'Type','line');
    %add callback fnc to the new line
    set(hline,'ButtonDownFcn',{@Line_Btn_down, hObject}); %do not pass handles -could be old
    handles.hSpec=hline;
end


if N>1
    cur = get(handles.SelectionBox,'Value');
    set(hline(cur),'Color',[1 0 0]);
end

% Update handles structure
guidata(hFig, handles);



%--------------------------------------------------------------------------
function update_plot2(hObject)
% PLE plot
hFig = ancestor(hObject,'figure','toplevel');
handles = guidata(hFig);
data = handles.data;
index = handles.filteredIndex;
file  = handles.selectedIndex;
yStep = handles.yStep;
if min(file)<=0
    return
end

mPx  = [];
xroi = [];
doCUT = 0;

if ~isempty(handles.markerPos)
    mPx =  handles.markerPos(1);
    doCUT = 1;
end

if ~isempty(handles.xroi)
    xroi = handles.xroi;
    doCUT = 2;
end



ifile=index(file);
N = length(ifile);
filename = handles.file_names(ifile);

%set CURRENT AXES to axes1 - main plot
set(handles.figure1,'CurrentAxes',handles.Axes2);
set(handles.Axes1,'NextPlot','replacechildren');

cutArr = get(handles.CutVariableDrop,'String');
cutStr = cutArr{get(handles.CutVariableDrop,'Value')};
rev = get(handles.PlotReversed,'Value');
nrm = get(handles.Normalize,   'Value');
bgd = get(handles.RemoveBgnd,  'Value');
tac = get(handles.ConsiderTac, 'Value');

if doCUT == 0
    return;
end


cut = zeros(N,2);
for i=1:N
    if rev
        j = N - i + 1;
    else
        j = i;
    end
    x = data(:,1,i);
    y = data(:,2,i);

    if tac 
        y= y/handles.tac(ifile(i));
    end
    
    if bgd
        miny = min(y); % alll yi and not: ..(intersect( find(xi>n(1)),find(xi<n(2)) ) ));
        y = y - miny;
        handles.bgd = miny;
    end

    %xPLE

    switch lower(cutStr)
        case 'power'
            if ~isempty(handles.pow)
                cut(i,1)= handles.pow(ifile(i));
            end
        case 'excitation'
            if ~isempty(handles.eng)
                cut(i,1)= handles.eng(ifile(i));
            end
            
        otherwise
            cut(i,1) = ifile(i); % like ''number'
    end

    %yPLE
    if (doCUT == 1 )                % use MARKER
        in1=find(x>mPx,1,'first');
        handles.cutEng = x(in1);
        handles.cutRng = [x(in1) x(in1)];
        cut(i,2) = y(in1);
    end
    if (doCUT == 2 )                % use XROI
        in1=find(x>min(xroi),1,'first');
        in2=find(x<max(xroi),1,'last');
        yr = y(in1:in2);
        handles.cutEng = (x(in1) +x(in2) )/2;
        handles.cutRng = [x(in1) x(in2)];
        cut(i,2) = sum(yr)./numel(yr);
    end
end
   
if (doCUT == 1 )
    handles.cutEng = x(in1);
    handles.cutRng = [x(in1) x(in1)];
elseif (doCUT == 2 )
    handles.cutEng = (x(in1) +x(in2) )/2;
    handles.cutRng = [x(in1) x(in2)];
end

plot(cut (:,1),cut(:,2),'-o');
% Update handles structure
guidata(hObject, handles);



%--------------------------------------------------------------------------
% --- Executes on MOUSE BTN press over line graph.
function Line_Btn_down(src, eventdata, hObject)
handles = guidata(hObject);
index = handles.filteredIndex;
file  = handles.selectedIndex;
line  = handles.selectedLine;
tool  = handles.tool;
cp = get(handles.Axes1,'CurrentPoint'); %read XY values

if file<=0
    return
end
doTool(tool, hObject, handles.Axes1);

% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function Axes1_ButtonDownFcn(hObject, eventdata, handles)
handles = guidata(hObject);
tool = handles.tool;

doTool(tool, hObject, handles.Axes1);

update_marker(hObject);
update_plot2(hObject);


%--------------------------------------------------------------------------
function doTool(tool, hObject,axes)
handles = guidata(hObject);
data = handles.data;
cp = get(axes,'CurrentPoint');

switch tool
    case 'zoom'
        %it will zoom by itself (built-in)
    case 'cut'
        if strcmp(get(handles.figure1,'SelectionType'),'alt')
            handles.markerPos=[];
        else
            handles.markerPos=cp(1,1:2);
            handles.xroi=[];
        end
    case 'cut ROI'
        % read doc rbbox - copied from there
        point1 = get(axes,'CurrentPoint');    % button down detected
        finalRect = rbbox;                   % return figure units
        point2 = get(axes,'CurrentPoint');    % button up detected
        point1 = point1(1,1:2);              % extract x and y
        point2 = point2(1,1:2);
        p1 = min(point1,point2);             % calculate locations
        offset = abs(point1-point2);         % and dimensions
        p2 = p1 + offset;

        %define Region of Interest
        xroi = [p1(1,1) p2(1,1)];
        %define dataR from xroi:
        xdata = data(:,1,1);                 % x-data of first spectrum
        if finalRect(1,3)~=0 && finalRect(1,4)~=0
            ileft  = find(xdata>xroi(1),1,'first');
            iright = find(xdata<xroi(2),1,'last');
            xr = [xdata(ileft) xdata(iright)];
            xr = xr(:);
            handles.xroi=xr;
            handles.markerPos=[];
        else
            handles.xroi=[];
            return
        end
    otherwise
        %
end
% Update handles structure
guidata(hObject, handles);





%% - SelectionBox -

% --- Executes during object creation, after setting all properties.
function SelectionBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SelectionBox_Callback(hObject, eventdata, handles)
hSpec = handles.hSpec;
cur = get(handles.SelectionBox,'Value');
set(hSpec,'Color',COLORPLOT);
set(hSpec(cur),'Color',[1 0 0]);



%% - YStepSlider -

% --- Executes during object creation, after setting all properties.
function YStepSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YStepSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function YStepSlider_Callback(hObject, eventdata, handles)
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

hFig = ancestor(hObject,'figure','toplevel');
val = get(hObject,'Value');
handles.yStep = val;

% Update handles structure
guidata(hFig, handles);
update_plot(hFig);



%% -PlotReversed -

% --- Executes on button press in PlotReversed.
function PlotReversed_Callback(hObject, eventdata, handles)
hFig = ancestor(hObject,'figure','toplevel');
val = get(hObject,'Value');
update_plot(hFig);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'KeyPressFcn',{@Figure_KeyPress,hObject,eventdata});


%--------------------------------------------------------------------------
% --- Executes on KBD key press when figure has focus.
function Figure_KeyPress(src, evnt, hObject, eventdata)

return


% --- Executes during object creation, after setting all properties.
function Axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Axes1�

%% - Tools Panel -
% --- Executes when selected object is changed in uiToolsPanel.
function uiToolsPanel_SelectionChangeFcn(source, eventdata, hObject)
handles = guidata(hObject);
sel = get(get(source,'SelectedObject'),'String');
handles.tool = sel;

hz = zoom(handles.figure1);
if strcmp(sel,'zoom')
    set(hz,'Enable','on');
else
    set(hz,'Enable','off');
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function uiToolsPanel_CreateFcn(hObject, eventdata, handles)
set(hObject,'SelectionChangeFcn',{@uiToolsPanel_SelectionChangeFcn,hObject});
none = findall(get(hObject,'Children'),'String','none');
set(hObject,'SelectedObject',none);


% --- Executes on button press in Normalize.
function Normalize_Callback(hObject, eventdata, handles)
hFig = ancestor(hObject,'figure','toplevel');
%handles = guidata(hFig);
val = get(hObject,'Value') ;
%guidata(hFig, handles);

set(handles.Axes1,'YLimMode','auto');

update_plot(hFig);
update_plot2(hFig);

% --- Executes on button press in RemoveBgnd.
function RemoveBgnd_Callback(hObject, eventdata, handles)
hFig = ancestor(hObject,'figure','toplevel');
val = get(hObject,'Value') ;
%guidata(hFig, handles);
update_plot(hFig);
update_plot2(hFig);

% --- Executes on button press in ConsiderTac.
function ConsiderTac_Callback(hObject, eventdata, handles)
hFig = ancestor(hObject,'figure','toplevel');
val = get(hObject,'Value');     % returns toggle state of ConsiderTac

if ~isempty(find(handles.tac==0,1) )
    %found at least one zero
    disp('PLE2: ERROR: at least one accumulation time is zero. ')
    set(hObject,'Value',0); 
end

update_plot(hFig);
update_plot2(hFig);



%% - Marker -

function update_marker(hObject)
%set CURRENT AXES to axes1 - main plot
%set(handles.figure1,'CurrentAxes',handles.axes1);
hFig = ancestor(hObject,'figure','toplevel');
handles = guidata(hFig);
hM = handles.hMarkers;
mP = handles.markerPos;
hS = handles.hSpec;
set(hFig,'CurrentAxes',handles.Axes1);
%marker must not extend beyond data or axes limits - get ranges
xdata = get(hS,'xdata');
ydata = get(hS,'ydata');
if iscell(xdata)
    xdata = cell2mat(xdata);
    ydata = cell2mat(ydata);
end
xm = [min(min(xdata)) max(max(xdata))];
ym = [min(min(ydata)) max(max(ydata))];
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
xdisp = [ max([xm(1) xlim(1)]) min([xm(2) xlim(2)]) ]; % max(lower values) min(higher values)
ydisp = [ max([ym(1) ylim(1)]) min([ym(2) ylim(2)]) ]; % max(lower values) min(higher values)
%set marker x,y
if ~isempty(handles.markerPos)
    tmp = get(gca,'ylim');
    set(hM(1,1),'Visible','on','Xdata',[mP(1,1) mP(1,1)],'Ydata',ydisp );
    tmp = get(gca,'xlim');
    set(hM(1,2),'Visible','on','Ydata',[mP(1,2) mP(1,2)],'Xdata',xdisp );
else
    set(hM(1,1),'Visible','off');
    set(hM(1,2),'Visible','off');
end
%handles.hMarkers = hM;
% Update handles structure
guidata(hObject, handles);



%% - Lower Buttons -



% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
dirname = fliplr(strtok(fliplr(pwd),filesep));
dirname = strcat(pwd,filesep,dirname,'_ple.mat');
[FileName,PathName,FilterIndex] = uiputfile(dirname) ;
if FileName == 0
    return
end
file_names  = handles.file_names;
eng = handles.eng;
pow = handles.pow;
tac = handles.tac;

save(FileName,'file_names','eng','pow','tac');


% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)

dirname = fliplr(strtok(fliplr(pwd),filesep));
dirname = strcat(pwd,filesep,dirname,'_ple.mat');
[FileName,PathName,FilterIndex] = uigetfile(dirname) ;
if FileName == 0
    return
end
load_matfile(FileName,hObject);


%--------------------------------------------------------------------------
%loads data from filename
function load_matfile(FileName, hObject)
% handles.selectedLine  = 1;      %selected line label (ie. X-) (+2 added in update_ bcs. header)
% handles.selectedIndex = 0;      %selected file in FileBox (+2 added in update_ to mark correct on bcs. headers)
% handles.filteredIndex = [];     %filtered index F- get real index R by: R=h.i_f(F)
handles.data = [];              %data of dsplayed files ( pts(2048) x 2 x num spec)
handles.eng  = [];
handles.pow  = [];
handles.tac  = [];
handles.tool = 'none';
% handles.hSpec= [];              % = handles to lines;
% handles.yStep = 0;              % step between spectra
% handles.hMarkers  = [];
% handles.markerPos = [];
handles.xroi = [];
handles.cutEng = 0;             % updated when PLE plot is drawn update_plot2()

load(FileName,'file_names','eng','pow','tac');
handles=guidata(hObject);
handles.file_names = [];
handles.field = [];
handles.energies = [];
handles.polarization = [];
handles.power = [];
check = 0;


%-----------
if exist('file_names','var')
    handles.file_names=file_names;
    handles.nFiles = size(file_names,1);
end
if exist('eng','var')
    handles.eng=eng;
else
    handles.eng=double(zeros(size(file_names,1)));
end
if exist('pow','var')
    handles.pow=pow;
else
    handles.pow=double(zeros(size(file_names,1)));
end
if exist('tac','var')
    handles.tac=tac;
else
    handles.tac=double(zeros(size(file_names,1)));
end

% Update handles structure
guidata(hObject, handles);
update_FileBox(hObject);





% --- Executes on button press in Export3D.
function Export3D_Callback(hObject, eventdata, handles)

hFig = ancestor(hObject,'figure','toplevel');
handles = guidata(hFig);
data  = handles.data;
index = handles.filteredIndex;
file  = handles.selectedIndex;

rev = get(handles.PlotReversed,'Value');
nrm = get(handles.Normalize,   'Value');
bgd = get(handles.RemoveBgnd,  'Value');

if min(file)<=0
    return
end
ifile=index(file);
N = length(ifile);

xlim = get(gca,'xlim');
xdata  = data(:,1,1);
ileft  = find(xdata>xlim(1),1,'first');
iright = find(xdata<xlim(2),1,'last');

X = data(ileft:iright,1,1) ; % less then 2048
Y = 1:N;

if bgd
    data = data - handles.bgd;
end

for i=1:N
    if nrm
        maxy = max( data(ileft:iright,2,i) );
        Z(:,i) = data(ileft:iright,2,i) ./ maxy;
    else
        Z(:,i) = data(ileft:iright,2,i);
    end
end

checkLabels = {'X','Y','Z'};
varNames = {'X','Y','Z'};
items = {X,Y,Z};
str = sprintf('Export to Workspace');
export2wsdlg(checkLabels,varNames,items,str);


% --- Executes on button press in MakePlot.
function MakePlot_Callback(hObject, eventdata, handles)
hFig = ancestor(hObject,'figure','toplevel');
handles = guidata(hFig);
index = handles.filteredIndex;
file  = handles.selectedIndex;

if min(file)<=0
    return
end
ifile=index(file);
powers = handles.pow(ifile);
energy = handles.eng(ifile);

ax1 = handles.Axes1;
nrm = get(handles.Normalize,   'Value');
cur = get(handles.SelectionBox,'Value');
   

fig = figure;
ax = copyobj(ax1,fig);
set(ax,'ButtonDownFcn',[]);
set(ax,'Units','normalized');
set(ax,'Position',[0.15 0.15 0.8 0.8]);
set(ax,'Box','on');

xlabel('Energy (meV)');
if nrm, ylabel('Normalised Intensity (arb.u.)');
else    ylabel('Intensity (arb.u.)');
end

cutArr = get(handles.CutVariableDrop,'String');
cutStr = cutArr{get(handles.CutVariableDrop,'Value')};
switch lower(cutStr)
    case 'power'
        str = strcat('P= ',mat2str(max(powers)),' - ',mat2str(min(powers)),' uW');
        annotation(fig, 'textbox', [0.2 0.8 0.2 0.1], 'String', str, ...
            'Margin', 0.001, 'FitHeightToText','on','LineStyle', 'none','Interpreter','none');
    case 'excitation'
        str = strcat('Lex= ',mat2str(max(energy)),' - ',mat2str(min(energy)),' nm');
        annotation(fig, 'textbox', [0.2 0.8 0.2 0.1], 'String', str, ...
            'Margin', 0.001, 'FitHeightToText', 'on', 'LineStyle', 'none', 'Interpreter', 'none');
    otherwise
        % like 'number'
end

if cur
    str=strcat(handles.file_names(ifile(cur)),' P= ',mat2str(handles.pow(ifile(cur))), ...
        'uW, Lex= ',mat2str(handles.eng(ifile(cur))),'nm');
    annotation(fig,'textbox', [0.2 0.7 0.5 0.1], 'String', str, ...
        'Color', [1 0 0], 'Margin', 0.001, 'FitHeightToText', 'on', 'LineStyle', 'none','Interpreter','none');
end

str = date;
annotation(fig,'textbox', [0.85 0.002 0.16 0.03], 'String',str, ...
    'Margin', 0.001, 'FitHeightToText','on', 'LineStyle', 'none');

str = pwd;
names = get(handles.SelectionBox,'String'); %list selected names
str = strcat(str,'/ ',' (', names(1,:), ' - ', names(end,:) ,') ');
annotation(fig,'textbox', [0.05 0.002 0.75 0.03], 'String',str, ...
    'Margin', 0.001, 'FitHeightToText','on', 'LineStyle', 'none','Interpreter','none');



% --- Executes on button press in MakePlotPLE.
function MakePlotPLE_Callback(hObject, eventdata, handles)

%%% %%%%%%%%%%%%%%%%% %%%
%%% TO DO : IMPLEMENT %%%
%%% %%%%%%%%%%%%%%%%% %%%



% --- Executes on button press in ExportPLE.
function ExportPLE_Callback(hObject, eventdata, handles)
hFig = ancestor(hObject,'figure','toplevel');
handles = guidata(hFig);
data = handles.data;
index = handles.filteredIndex;
file  = handles.selectedIndex;

ax = handles.Axes2;
hline=findobj(ax,'Type','line');
xdata = get(hline,'Xdata');
ydata = get(hline,'Ydata');

data = [xdata' ydata'];

str = sprintf('ple_%.3f', handles.cutEng);
str = strrep(str,'.','_');

cutArr = get(handles.CutVariableDrop,'String');
cutStr = cutArr{get(handles.CutVariableDrop,'Value')};

toexport = struct ('cut',handles.cutEng,'range',handles.cutRng,'xtype',cutStr,'x',xdata,'y',ydata);
checkLabels = {'PLE'};
varNames = {str};
items = {toexport};
str = sprintf('Export to Workspace');
export2wsdlg(checkLabels,varNames,items,str);



%--------------------------------------------------------------------------
% --- Executes on button press in LoadLog.
function LoadLog_Callback(hObject, eventdata, handles)
dirname = fliplr(strtok(fliplr(pwd),filesep));
dirname = strcat(pwd,filesep,dirname,'.log');
[FileName,PathName,FilterIndex] = uigetfile({'*.log','Log files';'*.*','All files'},'Load log',dirname);
if FileName==0
    return
end
fid = fopen(FileName);
num_files = size(handles.file_names,1);
tline = 0;

while tline ~= -1
    P   = [];
    Lex = [];
    tline = fgets(fid);

    name = strtok(tline,',');                           %get filename

    %get P value
    Ptxt=tline( strfind(tline,'P=') : size(tline,2));   %find 'P=' and cut all before
    Ptxt=strtok(Ptxt,',');                              %cut all after ','
    Pstr=sscanf(Ptxt,'P=%f%s');


    %get Lexc value
    Lex=tline( strfind(tline,'Lexc=') : size(tline,2)); %find 'Lexc=' and cut all before
    Lex=strtok(Lex,',');                                %cut all after ','
    Lex=sscanf(Lex,'Lexc=%fnm');

 %get t acc. value
    Tacc=tline( strfind(tline,'t=') : size(tline,2));   %find 't=' and cut all before
    Tacc=strtok(Tacc,',');                              %cut all after ','
    Tacc=sscanf(Tacc,'t=%fs');

    for j=1:num_files
        if strcmp(name,handles.file_names(j,1))         %look for filename
            if ~isempty(Pstr)
                P=Pstr(1);
                Pstr=char(Pstr(2:end))';
                switch Pstr
                    case 'uW'
                        P=P*1;
                    case 'nW'
                        P=P*0.001;
                    case 'mW'
                        P=P*1000;
                    otherwise
                        %
                end
                handles.pow(j,1)=P;                      %copy B to fields
            end

            if ~isempty(Lex)
                handles.eng(j,1) = Lex;
            end
            
            if ~isempty(Tacc)
                handles.tac(j,1) = Tacc;
            end

            break;
        end
    end
end
fclose(fid);
% Update handles structure
guidata(hObject, handles);
update_FileBox(hObject);




%% - Variable Drop -
% --- Executes on selection change in CutVariableDrop.
function CutVariableDrop_Callback(hObject, eventdata, handles)
hFig = ancestor(hObject,'figure','toplevel');
update_plot2(hFig);


% --- Executes during object creation, after setting all properties.
function CutVariableDrop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



