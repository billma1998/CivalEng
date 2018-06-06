function varargout = untitled(varargin)
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 04-Jun-2018 12:31:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
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


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function handles = getparametermeter(handles)
% Obtain the input to show the current Fo Number
handles.data.T_int = str2double(get(handles.T_int,'String'));
handles.data.T_top = str2double(get(handles.T_top,'String'));
handles.data.T_btm = str2double(get(handles.T_btm,'String'));
handles.data.T_lft = 25;
handles.data.T_rht = 25;
handles.data.L = str2double(get(handles.L,'String'));
handles.data.H = str2double(get(handles.H,'String'));
handles.data.dx = str2double(get(handles.dx,'String'));  
handles.data.dy = str2double(get(handles.dy,'String')); 
handles.data.tmax = str2double(get(handles.tmax,'String')); 
handles.data.dt = 0.01;
handles.data.epsilon = 0.001;

function conduction(handles)
% Numerically solve for transient conducton problem

% check if the Stop button is pressed: if not, proceed
set(handles.stop,'UserData',0);
% obtain the input parameters
L = handles.data.L;
H = handles.data.H;
dx = 0.01;
dy = 0.01;
tmax = 20;
dt =0.01;
epsilon = handles.data.epsilon;
r_x = handles.data.r_x;
r_y = handles.data.r_y;
% create the x, y meshgrid based on dx, dy
nx = uint32(L/dx + 1);
ny = uint32(H/dy + 1);
[X,Y] = meshgrid(linspace(0,L,nx),linspace(0,H,ny));
% take the center point of the domain
ic = uint32((nx-1)/2+1);
jc = uint32((ny-1)/2+1);   
% set initial and boundary conditions
T = handles.data.T_int*ones(ny,nx);
T(:,1) = handles.data.T_lft;
T(:,end) = handles.data.T_rht;
T(1,:) = handles.data.T_btm;
T(end,:) = handles.data.T_top;
Tmin = min(min(T));
Tmax = max(max(T));
% iteration, march in time
n = 0; 
nmax = uint32(tmax/dt);
while n < nmax
    if get(handles.stop,'UserData') == 1 
        break
    end
    n = n + 1;
    T_n = T;
    for j = 2:ny-1
        for i = 2:nx-1
            T(j,i) = T_n(j,i) + r_x*(T_n(j,i+1)-2*T_n(j,i)+T_n(j,i-1))...
                + r_y*(T_n(j+1,i)-2*T_n(j,i)+T_n(j-1,i));
        end
    end
    if uint16(n/50) == n/50 % refresh the plot every 50 time steps to save time     
        % plot temperature Tcontour
        handles.fig.cont = contourf(handles.Tcontour,X,Y,T,20);
        title(handles.Tcontour,sprintf('Time = %g s',n*dt)),
        colorbar('peer',handles.Tcontour),
        xlabel(handles.Tcontour,'x (m)'),ylabel(handles.Tcontour,'y (m)')
        axis(handles.Tcontour,'equal','tight'),
        % plot temperature at center point
        handles.fig.pl = scatter(handles.Tplot,n*dt,T(jc,ic),'r.');
        xlim(handles.Tplot,[0 tmax]),xlabel(handles.Tplot,'t (s)'),ylabel(handles.Tplot,'Tcenter')
        hold(handles.Tplot,'on')
        pause(0.01)
    end
    % check for convergence
    err = max(max(abs((T-T_n))));
    if err <= epsilon
        break
    end
   
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getparameter(handles);
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getparameter(handles);
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getparameter(handles);
% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getparameter(handles);
% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getparameter(handles);
% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
input = str2num(get(hObject,'String'));
if(isempty(input))
set(hObject,'String','0')
end
guidata(hObject,handles);
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getparametermeter(handles);

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
input = str2num(get(hObject,'String'));
if(isempty(input))
set(hObject,'String','0')
end
guidata(hObject,handles);
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.

function pushbutton1_Callback(hObject, eventdata, handles)

cla(handles.Tplot);
cla(handles.Tcontour,'reset');
handles = getparameter(handles);
conduction(handles);
guidata(hObject,handles);
% global a b c d e f
% A = get(handles.edit1,'String');
% B = get(handles.edit2,'String');
% C = get(handles.edit7,'String');
% D = get(handles.edit8,'String');
% E = get(handles.edit9,'String');
% F = get(handles.edit10,'String');
% a = str2num(A);
% b = str2num(B);
% c = str2num(C);
% d = str2num(D);
% e = str2num(E);
% f = str2num(F);
% clear;
% X = 500;
% PA = [0,0];
% PC = [(a-c)/2,b];
% PB = [0,500];
% PD = [(a-c)/2,X-b];
% PE = [0,X-b];
% PF = [0,X];
% PG = [a,X];
% PH = [a,X-b];
% PI = [(a+c)/2,X-b];
% PJ = [(a+c)/2,b];
% PK = [a,b];
% PL = [a,0];
% line(PA,PB)

%g = a + b+ c+d+e+f;
%h = num2str(g);
%set(handles.text8,'String',h);

% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
