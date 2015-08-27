function varargout = mock(varargin)
% MOCK MATLAB code for mock.fig
%      MOCK, by itself, creates a new MOCK or raises the existing
%      singleton*.
%
%      H = MOCK returns the handle to a new MOCK or the handle to
%      the existing singleton*.
%
%      MOCK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOCK.M with the given input arguments.
%
%      MOCK('Property','Value',...) creates a new MOCK or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mock_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mock_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mock

% Last Modified by GUIDE v2.5 08-Apr-2014 16:02:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mock_OpeningFcn, ...
                   'gui_OutputFcn',  @mock_OutputFcn, ...
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

% --- Executes just before mock is made visible.
function mock_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mock (see VARARGIN)

% Choose default command line output for mock
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes mock wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mock_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function TnDensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TnDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TnDensity_Callback(hObject, eventdata, handles)
% hObject    handle to TnDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TnDensity as text
%        str2double(get(hObject,'String')) returns contents of TnDensity as a double
TnDensity = str2double(get(hObject, 'String'));
if isnan(TnDensity)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new TnDensity value
handles.metricdata.TnDensity = TnDensity;
guidata(hObject,handles)

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

handles.metricdata.TnDensity = 1;
handles.metricdata.XBDensity = 1;
handles.metricdata.pCa = 4;
handles.metricdata.StartLength = 1150;
handles.metricdata.kxscaler = 3;
handles.metricdata.L_TITIN = 247;

handles.metricdata.TnKOType = 0;
handles.metricdata.XBKOType = 0;

set(handles.TnDensity, 'String', handles.metricdata.TnDensity);
set(handles.XBDensity,  'String', handles.metricdata.XBDensity);
set(handles.pCaEdit, 'String', handles.metricdata.pCa);
set(handles.StartLengthEdit,  'String', handles.metricdata.StartLength);
set(handles.kxscalerEdit, 'String', handles.metricdata.kxscaler);
set(handles.LTitinEdit,  'String', handles.metricdata.L_TITIN);

set(handles.TnRandom, 'Value', 1 - handles.metricdata.TnKOType);
set(handles.TnUniform, 'Value', handles.metricdata.TnKOType);
set(handles.XBRandom, 'Value', 1 - handles.metricdata.XBKOType);
set(handles.XBUniform, 'Value', handles.metricdata.XBKOType);

% Update handles structure
guidata(handles.figure1, handles);


% --- Executes on button press in XBUniform.
function XBUniform_Callback(hObject, eventdata, handles)
% hObject    handle to XBUniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of XBUniform

set(handles.XBRandom, 'Value', 1 - get(hObject,'Value'));
handles.metricdata.XBKOType = 1;
guidata(hObject, handles)


% --- Executes on button press in XBRandom.
function XBRandom_Callback(hObject, eventdata, handles)
% hObject    handle to XBRandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of XBRandom

set(handles.XBUniform, 'Value', 1 - get(hObject,'Value'));
handles.metricdata.XBKOType = 0;
guidata(hObject, handles)


function XBDensity_Callback(hObject, eventdata, handles)
% hObject    handle to XBDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XBDensity as text
%        str2double(get(hObject,'String')) returns contents of XBDensity as a double
XBDensity = str2double(get(hObject, 'String'));
if isnan(XBDensity)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new TnDensity value
handles.metricdata.XBDensity = XBDensity;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function XBDensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XBDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TnRandom.
function TnRandom_Callback(hObject, eventdata, handles)
% hObject    handle to TnRandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TnRandom
set(handles.TnUniform, 'Value', 1 - get(hObject,'Value'));
handles.metricdata.TnKOType = 0;
guidata(hObject, handles)


% --- Executes on button press in TnUniform.
function TnUniform_Callback(hObject, eventdata, handles)
% hObject    handle to TnUniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TnUniform
set(handles.TnRandom, 'Value', 1 - get(hObject,'Value'));
handles.metricdata.TnKOType = 1;
guidata(hObject, handles)


function pCaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pCaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pCaEdit as text
%        str2double(get(hObject,'String')) returns contents of pCaEdit as a double
pCa = str2double(get(hObject, 'String'));
if isnan(pCa)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new TnDensity value
handles.metricdata.pCa = pCa;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function pCaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pCaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kxscalerEdit_Callback(hObject, eventdata, handles)
% hObject    handle to kxscalerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kxscalerEdit as text
%        str2double(get(hObject,'String')) returns contents of kxscalerEdit as a double
kxscaler = str2double(get(hObject, 'String'));
if isnan(kxscaler)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new TnDensity value
handles.metricdata.kxscaler = kxscaler;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function kxscalerEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kxscalerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StartLengthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to StartLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartLengthEdit as text
%        str2double(get(hObject,'String')) returns contents of StartLengthEdit as a double
StartLength = str2double(get(hObject, 'String'));
if isnan(StartLength)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new TnDensity value
handles.metricdata.StartLength = StartLength;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function StartLengthEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartLengthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LTitinEdit_Callback(hObject, eventdata, handles)
% hObject    handle to LTitinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LTitinEdit as text
%        str2double(get(hObject,'String')) returns contents of LTitinEdit as a double
L_TITIN = str2double(get(hObject, 'String'));
if isnan(L_TITIN)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new TnDensity value
handles.metricdata.L_TITIN = L_TITIN;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function LTitinEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LTitinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GoButton.
function GoButton_Callback(hObject, eventdata, handles)
% hObject    handle to GoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
init_params;

RequestedRuns = 4;
filaments.L_TITIN = handles.metricdata.L_TITIN;
StartLength = handles.metricdata.StartLength;
StiffScale.kxscaler = handles.metricdata.kxscaler;
pCa = handles.metricdata.pCa;

knockout.TnKOType = handles.metricdata.TnKOType;
knockout.XBKOType = handles.metricdata.XBKOType;
knockout.TnDensity = handles.metricdata.TnDensity;
knockout.XBDensity = handles.metricdata.XBDensity;

[Steps, Means, Vars, IndexThalf, Binder] = RunSeveral(RequestedRuns, DataParams, StartLength, pCa, StiffScale, filaments, knockout, coop, TFRateScale, tcparam);

clf(figure(1))
SimpleForceTime(DataParams.dt, Steps(1,:), Steps(2,:))

disp('done')
