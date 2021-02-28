function varargout = rearrange(varargin)
% REARRANGE M-file for rearrange.fig
% MB 18/05/2006
% ver 1.0
%      REARRANGE, by itself, creates a new REARRANGE or raises the existing
%      singleton*.
%
%      H = REARRANGE returns the handle to a new REARRANGE or the handle to
%      the existing singleton*.
%
%      REARRANGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REARRANGE.M with the given input arguments.
%
%      REARRANGE('Property','Value',...) creates a new REARRANGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rearrange_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rearrange_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rearrange

% Last Modified by GUIDE v2.5 18-May-2006 12:39:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rearrange_OpeningFcn, ...
                   'gui_OutputFcn',  @rearrange_OutputFcn, ...
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


% --- Executes just before rearrange is made visible.
function rearrange_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rearrange (see VARARGIN)

% Choose default command line output for rearrange
handles.output = hObject;
handles.list={};
for k=1:size(varargin{1},2)
    handles.list{k}=char(varargin{1}(k));
end
% Update handles structure
guidata(hObject, handles);
%update listbox:
set(handles.listbox1,'String',handles.list);

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
% FigPos=get(0,'DefaultFigurePosition');
% FigWidth=215;FigHeight=88;
% if isempty(gcbf)
%     ScreenUnits=get(0,'Units');
%     set(0,'Units','points');
%     ScreenSize=get(0,'ScreenSize');
%     set(0,'Units',ScreenUnits);
% 
%     FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
%     FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
% else
%     GCBFOldUnits = get(gcbf,'Units');
%     set(gcbf,'Units','points');
%     GCBFPos = get(gcbf,'Position');
%     set(gcbf,'Units',GCBFOldUnits);
%     FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
%                    (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
% end
% FigPos(3:4)=[FigWidth FigHeight];
% set(hObject, 'position', FigPos);

% UIWAIT makes rearrange wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = rearrange_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
 %end
 
 
 
% --- Executes on button press in bt_up.
function bt_up_Callback(hObject, eventdata, handles)
% hObject    handle to bt_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imax=size(handles.list,2);
index=get(handles.listbox1,'Value');
if isequal(index,1)
    % do nothing
else
    %switch with index-1
    str=handles.list(1,index);
    handles.list(:,index)=handles.list(:,index-1);
    handles.list(:,index-1)=str;
    %update listbox:
    set(handles.listbox1,'String',handles.list);
    set(handles.listbox1,'Value',index-1);
    guidata(hObject, handles);
end



% --- Executes on button press in bt_dn.
function bt_dn_Callback(hObject, eventdata, handles)
% hObject    handle to bt_dn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imax=size(handles.list,2);
index=get(handles.listbox1,'Value');
if isequal(index,imax)
    % do nothing
else
    %switch with index-1
    str=handles.list(1,index);
    handles.list(:,index)=handles.list(:,index+1);
    handles.list(:,index+1)=str;
    %update listbox:
    set(handles.listbox1,'String',handles.list);
    set(handles.listbox1,'Value',index+1);
    guidata(hObject, handles);
end





% --- Executes on button press in bt_cancel.
function bt_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to bt_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output='';
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in bt_OK.
function bt_OK_Callback(hObject, eventdata, handles)
% hObject    handle to bt_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output=handles.list;
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end







% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


