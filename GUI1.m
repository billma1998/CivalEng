function varargout = GUI1(varargin)
% GUI1 MATLAB code for GUI1.fig
%      GUI1, by itself, creates a new GUI1 or raises the existing
%      singleton*.
%
%      H = GUI1 returns the handle to a new GUI1 or the handle to
%      the existing singleton*.
%
%      GUI1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI1.M with the given input arguments.
%
%      GUI1('Property','Value',...) creates a new GUI1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI1

% Last Modified by GUIDE v2.5 07-Jun-2018 11:40:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI1_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI1_OutputFcn, ...
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

function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = handles.sel_data;
switch s
    case 1
      func1
    case 2
      func2
    case 3
     func3   
    case 4
   func4
    case 5
   func5 
    case 6
      func6
end
function diergeshu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diergeshu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a b c d e
X = get(handles.diyigeshu,'String');
Y = get(handles.diergeshu,'String');

P  = get(handles.thermaldelta,'String');
G = get(handles.densitydelta,'String');
H = get(handles.heatdelta,'String');
a = str2num(X);
b = str2num(Y);
c = str2num(P); 
d = str2double(G);
e = str2double(H);

function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = {'cross section','beam','column','frame','animation of Temp in cross section','6'};
[sel,ok] = listdlg('ListString',str,'Name','Choose the function','PromptString','Select','SelectionMode','single');
handles.sel_data = sel
guidata(hObject,handles)




function func1(handles)
    global a b c d e
%     x = get(handles.diyigeshu,'String');
%     y = get(handles.diergeshu,'String');
%     a = str2double(x);
%    
%     b = str2double(y);
%     c = a+b
%     d = num2str(c);
% %     set(handles.jieguo,'String',d);
%     guidata(hObject.handles);
%  a = str2double(get(handles.diyigeshu,'String'));
% b =str2double(get(handles.diergeshu,'String'));
A = createpde('thermal','steadystate');
% point1 = [3 4 -.075 .075 .075 -.075  -.1125 -.1125 .1125 .1125];
% point2 = [3 4 -.075 -.0075 -.0075 -.075  -.0975 -.0975 .0975 .0975];
% point3 = [3 4 .075 .0075 .0075 .075  -.0975 -.0975 .0975 .0975];

 point1 = [3 4 -.075 .075 .075 -.075  -.1125 -.1125 .1125 .1125];
 point2 = [3 4 -.075 -.0075 -.0075 -.075  -.0975 -.0975 .0975 .0975];
 point3 = [3 4 .075 .0075 .0075 .075  -.0975 -.0975 .0975 .0975];

reverse = [point1; point2;point3]';
g = decsg(reverse,'R1-R2-R3',['R1'; 'R2';'R3']');
geometryFromEdges(A,g);
 
 
thermalBC(A,'Edge',2,'HeatFlux',-22.2);
thermalBC(A,'Edge',3,'HeatFlux',-33.3);
thermalBC(A,'Edge',4,'HeatFlux',-33.3);
thermalBC(A,'Edge',5,'HeatFlux',-2.22);
thermalBC(A,'Edge',6,'HeatFlux',-2.22);
thermalBC(A,'Edge',7,'HeatFlux',-21.09);
thermalBC(A,'Edge',8,'HeatFlux',-21.09);
thermalBC(A,'Edge',9,'HeatFlux',-21.09);
thermalBC(A,'Edge',10,'HeatFlux',-21.09);
thermalBC(A,'Edge',11,'HeatFlux',-2.22);
thermalBC(A,'Edge',12,'HeatFlux',-2.22);
thermalBC(A,'Edge',1,'Temperature',a);
 
thermalProperties(A,'ThermalConductivity',e);
 
 
generateMesh(A,'Hmax',0.2);
 
R = solve(A);
T = R.Temperature;
figure
pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet');
axis equal
title 'Temperature, Steady State Solution'
 
 
B = createpde('thermal','transient');
 
point1 = [3 4 -.075 .075 .075 -.075  -.1125 -.1125 .1125 .1125];
point2 = [3 4 -.075 -.0075 -.0075 -.075  -.0975 -.0975 .0975 .0975];
point3 = [3 4 .075 .0075 .0075 .075  -.0975 -.0975 .0975 .0975];
gdm = [point1; point2;point3]';
g = decsg(gdm,'R1-R2-R3',['R1'; 'R2';'R3']');
geometryFromEdges(B,g);
 
 
thermalProperties(B,'MassDensity', c,'SpecificHeat',d, 'ThermalConductivity',e);
                             
                                
  
thermalBC(B,'Edge',2,'HeatFlux',-22.2);
thermalBC(B,'Edge',3,'HeatFlux',-33.3);
thermalBC(B,'Edge',4,'HeatFlux',-33.3);
thermalBC(B,'Edge',5,'HeatFlux',-2.22);
thermalBC(B,'Edge',6,'HeatFlux',-2.22);
thermalBC(B,'Edge',7,'HeatFlux',-21.09);
thermalBC(B,'Edge',8,'HeatFlux',-21.09);
thermalBC(B,'Edge',9,'HeatFlux',-21.09);
thermalBC(B,'Edge',10,'HeatFlux',-21.09);
thermalBC(B,'Edge',11,'HeatFlux',-2.22);
thermalBC(B,'Edge',12,'HeatFlux',-2.22);
thermalBC(B,'Edge',1,'Temperature',a);
 
msh = generateMesh(B,'Hmax',0.2);
 
tlist = 0:.1:50;
thermalIC(B,b);
R = solve(B,tlist);
T = R.Temperature;
 
 
getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
 
[~,nid] = getClosestNode( msh.Nodes, .5, 0 );
 
 
h = figure;
h.Position = [1 1 2 1].*h.Position;
subplot(1,2,1);
axis equal
pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet');
axis equal
title 'transient solution'
subplot(1,2,2);
axis equal
plot(tlist, T(nid,:));
grid on
title 'Temperature as different Time C';
xlabel 'seconds'
ylabel 'Temperature C'
% close all
% F(length(tlist)) = struct('cdata',[],'colormap',[]);
% for tt=1:length(tlist)
% pdeplot(A,'XYData',T(:,tt),'Contour','on','ColorMap','jet');
% drawnow;
% F(tt) = getframe;
% end

 
 
 
% k = @(~,state) 0.3+0.003*state.u;
%  
% thermalProperties(A,'ThermalConductivity',k);
%  
% R = solve(A);
% T = R.Temperature;
% figure
% pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet'); 
% axis equal
% title 'Temperature, Steady State Solution'
%  
%  
% thermalProperties(B,'ThermalConductivity',k,'MassDensity',1,'SpecificHeat',1);
%  
% thermalIC(B,0);
% R = solve(B,tlist);
% T = R.Temperature;
%  
%  
% h = figure;
% h.Position = [1 1 2 1].*h.Position;
% subplot(1,2,1);
% axis equal
% pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet');
% axis equal
% title 'Temperature, Final Time, Transient Solution'
% subplot(1,2,2);
% axis equal
% plot(tlist(1:size(T,2)), T(nid,:));
% grid on
% title 'Temperature at Right Edge as a Function of Time (Nonlinear)';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'
%  
 
%  
% %for animation
% close all
% F(length(tlist)) = struct('cdata',[],'colormap',[]);
% for tt=1:length(tlist)
% pdeplot(A,'XYData',T(:,tt),'Contour','on','ColorMap','jet');
% drawnow;
% F(tt) = getframe;
% 
% end


function func2(handles)
global a b c d e
A = createpde('thermal');
g = decsg([3 4 -1.5 1.5 1.5 -1.5 0 0 .2 .2]');
geometryFromEdges(A,g);

con = c; 
den = d; 
capa = e; 
supp = 500*a; 
FCON = @(region,state) con*region.y;
FCAPA = @(region,state) capa*region.y;
FSUPP = @(region,state) supp*region.y;

thermalProperties(A,'ThermalConductivity',FCON);
internalHeatSource(A,FSUPP);


thermalBC(A,'Edge',2,'Temperature',100);

coff = @(region,~) 50*region.y;
thermalBC(A,'Edge',3,'ConvectionCoefficient',coff,'AmbientTemperature',100);
                   
leftHF = @(region,~) 5000*region.y;
thermalBC(A,'Edge',4,'HeatFlux',leftHF);


generateMesh(A,'Hmax',0.2);

result = solve(A);
T = result.Temperature;
figure; 
pdeplot(A,'XYData',T,'Contour','on'); 
axis equal
title 'Steady State Temperature';

B = createpde('thermal','transient');

g = decsg([3 4 -1.5 1.5 1.5 -1.5 0 0 .2 .2]');
geometryFromEdges(B,g);

generateMesh(B,'Hmax',0.2);

thermalProperties('SpecificHeat',FCAPA ,B,'ThermalConductivity',FCON,'MassDensity',den);                           

internalHeatSource(B,FSUPP);

thermalBC(B,'Edge',2,'Temperature',100);
thermalBC(B,'Edge',3,'ConvectionCoefficient',coff,'AmbientTemperature',b);
thermalBC(B,'Edge',4,'HeatFlux',leftHF);

tfinal = 20000;
tlist = 0:100:tfinal;
thermalIC(B,b);
B.SolverOptions.ReportStatistics = 'on';

result = solve(B,tlist);
T = result.Temperature;

figure; 
pdeplot(B,'XYData',T(:,end),'Contour','on'); 
axis equal
title(sprintf('Transient Temperature at Final Time (%g seconds)',tfinal));

p = B.Mesh.Nodes;
nodesLeftEnd  = find(p(1,:) < -1.5+eps);
nodeCenter = nodesLeftEnd(p(2,nodesLeftEnd) < eps);
nodeOuter = nodesLeftEnd(p(2,nodesLeftEnd) > 0.2-eps);

figure; 
plot(tlist,T(nodeCenter,:)); 
hold all 
plot(tlist,T(nodeOuter,:),'--');
title 'Temperature at Left End as a Function of Time'
xlabel 'Time, seconds'
ylabel 'Temperature C'
grid on;
legend('Center Axis','Outer Surface','Location','SouthEast');


function func3(handles);
global a b c d e 
A = createpde('thermal','steadystate');
g = decsg([3 4 -0.5 0.5 0.5 -0.5 -4.5 -4.5 4.5 4.5]');
geometryFromEdges(A,g);
 


 
thermalBC(A,'Edge',2,'HeatFlux',-22.2);
thermalBC(A,'Edge',3,'HeatFlux',-33.3);
thermalBC(A,'Edge',4,'HeatFlux',-33.3);
% thermalBC(A,'Edge',5,'HeatFlux',-2.22);
% thermalBC(A,'Edge',6,'HeatFlux',-2.22);
% thermalBC(A,'Edge',7,'HeatFlux',-21.09);
% thermalBC(A,'Edge',8,'HeatFlux',-21.09);
% thermalBC(A,'Edge',9,'HeatFlux',-21.09);
% thermalBC(A,'Edge',10,'HeatFlux',-21.09);
% thermalBC(A,'Edge',11,'HeatFlux',-2.22);
% thermalBC(A,'Edge',12,'HeatFlux',-2.22);
thermalBC(A,'Edge',1,'Temperature',a);
 
thermalProperties(A,'ThermalConductivity',e);
 
 
generateMesh(A,'Hmax',0.2);
 
R = solve(A);
T = R.Temperature;
figure
pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet');
axis equal
title 'Temperature, Steady State Solution'
 
 
B = createpde('thermal','transient');

g = decsg([3 4 -0.5 0.5 0.5 -0.5 -4.5 -4.5 4.5 4.5]');

geometryFromEdges(B,g);
% pdegplot(B,'EdgeLabels','on');
% axis([-9 9 -9 9]);


 
thermalProperties(B,'MassDensity',d, 'ThermalConductivity',c,'SpecificHeat',e);
                             
                                
 
 
thermalBC(B,'Edge',2,'HeatFlux',-10);
thermalBC(B,'Edge',1,'Temperature',a);
 
msh = generateMesh(B,'Hmax',0.02);
 
tlist = 0:10000:300000;
thermalIC(B,b);
R = solve(B,tlist);
T = R.Temperature;
 
 
getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
 
[~,nid] = getClosestNode( msh.Nodes, .5, 0 );
h = figure;
h.Position = [1 1 2 1].*h.Position;
subplot(1,2,1); 
axis equal  
pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet'); 
axis equal
title 'Temperature, Final Time, Transient Solution'
subplot(1,2,2); 
axis equal
plot(tlist, T(nid,:)); 
grid on
title 'Temperature at Right Edge as a Function of Time';
xlabel 'Time, seconds'
ylabel 'Temperature, degrees-Celsius'

 
 

 

function diyigeshu_Callback(hObject, eventdata, handles)
% hObject    handle to diyigeshu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
input = str2num(get(hObject,'String'));
if(isempty(input))
    set(hObject,'String',0)
end
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of diyigeshu as text
%        str2double(get(hObject,'String')) returns contents of diyigeshu as a double


% --- Executes during object creation, after setting all properties.
function diyigeshu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diyigeshu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diergeshu_Callback(hObject, eventdata, handles)
% hObject    handle to diergeshu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
input = str2num(get(hObject,'String'));
if(isempty(input))
    set(hObject,'String',0)
end
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of diergeshu as text
%        str2double(get(hObject,'String')) returns contents of diergeshu as a double


% --- Executes during object creation, after setting all properties.





function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = questdlg('close this windows?','MATLAB quest')
if answer == 'Yes'
    close
end

 
 

function func4(handles)
global a b c d e

A = createpde('thermal','steadystate');
point2 = [2 3 -8 0 8 1 9 1];
point1 = [2 3 -10 0 10 0 10 0];

gdm = [point1; point2]';
g = decsg(gdm,'R1-R2',['R1'; 'R2']');
geometryFromEdges(A,g);


thermalBC(A,'Edge',2,'HeatFlux',-10);
thermalBC(A,'Edge',1,'Temperature',a);

thermalProperties(A,'ThermalConductivity',c);


generateMesh(A,'Hmax',0.2);

R = solve(A);
T = R.Temperature;
figure
pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet'); 
axis equal
title 'Temperature, Steady State Solution'
B = createpde('thermal','transient');

r1 = [2 3 -10 0 10 0 10 0];
r2 = [2 3 -8 0 8 1 9 1];
gdm = [r1; r2]';
g = decsg(gdm,'R1-R2',['R1'; 'R2']');
geometryFromEdges(B,g);

thermalProperties(B,'ThermalConductivity',c,'MassDensity',d,'SpecificHeat',e);

thermalBC(B,'Edge',1,'HeatFlux',-10);
% thermalBC(B,'Edge',6,'Temperature',@transientBCHeatedBlock);
thermalBC(B,'Edge',1,'Temperature',a);
msh = generateMesh(B,'Hmax',0.2);
tlist = 0:1000:400000;
thermalIC(B,b);
R = solve(B,tlist);
T = R.Temperature;


getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);

[~,nid] = getClosestNode( msh.Nodes, .5, 0 );


h = figure;
h.Position = [1 1 2 1].*h.Position;
subplot(1,2,1); 
axis equal  
pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet'); 
axis equal
title 'Temperature, Final Time, Transient Solution'
subplot(1,2,2); 
axis equal
plot(tlist, T(nid,:)); 
grid on
title 'Temperature at Right Edge as a Function of Time';
xlabel 'Time, seconds'
ylabel 'Temperature, degrees-Celsius'
% A = createpde('thermal','steadystate');
% point2 = [2 3 -8 0 8 1 9 1];
% point1 = [2 3 -10 0 10 0 10 0];
% 
% gdm = [point1; point2]';
% g = decsg(gdm,'R1-R2',['R1'; 'R2']');
% geometryFromEdges(A,g);
% 
% 
% thermalBC(A,'Edge',2,'HeatFlux',-10);
% thermalBC(A,'Edge',1,'Temperature',a);
% 
% thermalProperties(A,'ThermalConductivity',51.5);
% 
% 
% generateMesh(A,'Hmax',0.2);
% 
% R = solve(A);
% T = R.Temperature;
% figure
% pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet'); 
% axis equal
% title 'Temperature, Steady State Solution'
% 
% 
% B = createpde('thermal','transient');
% 
% point2 = [2 3 -8 0 8 1 9 1];
% point1 = [2 3 -10 0 10 0 10 0];
% 
% gdm = [point1; point2]';
% g = decsg(gdm,'R1-R2',['R1'; 'R2']');
% geometryFromEdges(B,g);
% 
% 
% thermalProperties(B,'MassDensity',7800, 'ThermalConductivity',51.5,'SpecificHeat',483);
%                               
%                                 
% thermalBC(B,'Edge',2,'HeatFlux',-10);
% thermalBC(B,'Edge',1,'Temperature',@transientBCHeatedBlock);
% 
% msh = generateMesh(B,'Hmax',0.2);
% 
% tlist = 0:.1:10;
% thermalIC(B,0);
% R = solve(B,tlist);
% T = R.Temperature;
% 
% 
% getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
% 
% [~,nid] = getClosestNode( msh.Nodes, .5, 0 );
% 
% 
% h = figure;
% h.Position = [1 1 2 1].*h.Position;
% subplot(1,2,1); 
% axis equal
% pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet'); 
% axis equal
% title 'Temperature, Final Time, Transient Solution'
% subplot(1,2,2); 
% axis equal
% plot(tlist, T(nid,:)); 
% grid on
% title 'Temperature at Right Edge as a Function of Time';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'



% k = @(~,state) 0.3+0.003*state.u;
% 
% thermalProperties(A,'ThermalConductivity',k);
%  
% 
% R = solve(A);
% T = R.Temperature;
% figure
% pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet');  
% axis equal
% title 'Temperature, Steady State Solution'
% 
% 
% thermalProperties(B,'ThermalConductivity',k,'MassDensity',7800,'SpecificHeat',482.5);
%  
% 
% thermalIC(B,b);
% R = solve(B,tlist);
% T = R.Temperature;
% 
% 
% h = figure;
% h.Position = [1 1 2 1].*h.Position;
% subplot(1,2,1); 
% axis equal
% pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet'); 
% axis equal
% title 'Temperature, Final Time, Transient Solution'
% subplot(1,2,2); 
% axis equal
% plot(tlist(1:size(T,2)), T(nid,:)); 
% grid on
% title 'Temperature at Right Edge as a Function of Time (Nonlinear)';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'

% --- Executes just before GUI1 is made visible.
function GUI1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI1 (see VARARGIN)

% Choose default command line output for GUI1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% function varargout = GUI1(varargin)
% GUI1 MATLAB code for GUI1.fig
%      GUI1, by itself, creates a new GUI1 or raises the existing
%      singleton*.
% 
%      H = GUI1 returns the handle to a new GUI1 or the handle to
%      the existing singleton*.
% 
%      GUI1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI1.M with the given input arguments.
% 
%      GUI1('Property','Value',...) creates a new GUI1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI1_OpeningFcn via varargin.
% 
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
% 
% See also: GUIDE, GUIDATA, GUIHANDLES
% 
% Edit the above text to modify the response to help GUI1
% 
% Last Modified by GUIDE v2.5 05-Jun-2018 23:49:50
% 
% Begin initialization code - DO NOT EDIT
% gui_Singleton = 1;
% gui_State = struct('gui_Name',       mfilename, ...
%                    'gui_Singleton',  gui_Singleton, ...
%                    'gui_OpeningFcn', @GUI1_OpeningFcn, ...
%                    'gui_OutputFcn',  @GUI1_OutputFcn, ...
%                    'gui_LayoutFcn',  [] , ...
%                    'gui_Callback',   []);
% if nargin && ischar(varargin{1})
%     gui_State.gui_Callback = str2func(varargin{1});
% end
% 
% if nargout
%     [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
% else
%     gui_mainfcn(gui_State, varargin{:});
% end
% End initialization code - DO NOT EDIT
% 
% function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% s = handles.sel_data;
% switch s
%     case 1
%     func1
%     case 2
%       func2
%     case 3
%      func3   
%     case 4
%    func4
% end
% 
% function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% str = {'cross section','beam','column','frame'};
% [sel,ok] = listdlg('ListString',str,'Name','Choose the function','PromptString','Select','SelectionMode','single');
% handles.sel_data = sel
% guidata(hObject,handles)
% 
% 
% function func1(handles)
%     clear all
% A = createpde('thermal','steadystate');
% point1 = [3 4 -.075 .075 .075 -.075  -.1125 -.1125 .1125 .1125];
% point2 = [3 4 -.075 -.0075 -.0075 -.075  -.0975 -.0975 .0975 .0975];
% point3 = [3 4 .075 .0075 .0075 .075  -.0975 -.0975 .0975 .0975];
% 
%  point1 = [3 4 -.075 .075 .075 -.075  -.1125 -.1125 .1125 .1125];
%  point2 = [3 4 -.075 -.0075 -.0075 -.075  -.0975 -.0975 .0975 .0975];
%  point3 = [3 4 .075 .0075 .0075 .075  -.0975 -.0975 .0975 .0975];
% 
% reverse = [point1; point2;point3]';
% g = decsg(reverse,'R1-R2-R3',['R1'; 'R2';'R3']');
% geometryFromEdges(A,g);
%  
%  
% thermalBC(A,'Edge',2,'HeatFlux',-22.2);
% thermalBC(A,'Edge',3,'HeatFlux',-33.3);
% thermalBC(A,'Edge',4,'HeatFlux',-33.3);
% thermalBC(A,'Edge',5,'HeatFlux',-2.22);
% thermalBC(A,'Edge',6,'HeatFlux',-2.22);
% thermalBC(A,'Edge',7,'HeatFlux',-21.09);
% thermalBC(A,'Edge',8,'HeatFlux',-21.09);
% thermalBC(A,'Edge',9,'HeatFlux',-21.09);
% thermalBC(A,'Edge',10,'HeatFlux',-21.09);
% thermalBC(A,'Edge',11,'HeatFlux',-2.22);
% thermalBC(A,'Edge',12,'HeatFlux',-2.22);
% thermalBC(A,'Edge',1,'Temperature',400);
%  
% thermalProperties(A,'ThermalConductivity',51.5);
%  
%  
% generateMesh(A,'Hmax',0.2);
%  
% R = solve(A);
% T = R.Temperature;
% figure
% pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet');
% axis equal
% title 'Temperature, Steady State Solution'
%  
%  
% B = createpde('thermal','transient');
%  
% point1 = [3 4 -.075 .075 .075 -.075  -.1125 -.1125 .1125 .1125];
% point2 = [3 4 -.075 -.0075 -.0075 -.075  -.0975 -.0975 .0975 .0975];
% point3 = [3 4 .075 .0075 .0075 .075  -.0975 -.0975 .0975 .0975];
% gdm = [point1; point2;point3]';
% g = decsg(gdm,'R1-R2-R3',['R1'; 'R2';'R3']');
% geometryFromEdges(B,g);
%  
%  
% thermalProperties(B,'MassDensity',7800, 'ThermalConductivity',51.5,'SpecificHeat',482.5);
%                              
%                                 
%  
% thermalBC(B,'Edge',2,'AmbientTemperature', b); 
% thermalBC(B,'Edge',2,'HeatFlux',-22.2);
% thermalBC(B,'Edge',3,'HeatFlux',-33.3);
% thermalBC(B,'Edge',4,'HeatFlux',-33.3);
% thermalBC(B,'Edge',5,'HeatFlux',-2.22);
% thermalBC(B,'Edge',6,'HeatFlux',-2.22);
% thermalBC(B,'Edge',7,'HeatFlux',-21.09);
% thermalBC(B,'Edge',8,'HeatFlux',-21.09);
% thermalBC(B,'Edge',9,'HeatFlux',-21.09);
% thermalBC(B,'Edge',10,'HeatFlux',-21.09);
% thermalBC(B,'Edge',11,'HeatFlux',-2.22);
% thermalBC(B,'Edge',12,'HeatFlux',-2.22);
% thermalBC(B,'Edge',1,'Temperature',400);
%  
% msh = generateMesh(B,'Hmax',0.2);
%  
% tlist = 0:.1:10;
% thermalIC(B,b);
% R = solve(B,tlist);
% T = R.Temperature;
%  
%  
% getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
%  
% [~,nid] = getClosestNode( msh.Nodes, .5, 0 );
%  
%  
% h = figure;
% h.Position = [1 1 2 1].*h.Position;
% subplot(1,2,1);
% axis equal
% pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet');
% axis equal
% title 'Temperature, Final Time, Transient Solution'
% subplot(1,2,2);
% axis equal
% plot(tlist, T(nid,:));
% grid on
% title 'Temperature at Right Edge as a Function of Time';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'

function AppHelp(varargin)
% Create help message
dlgname = 'About GUI';
txt = {'First input box is for the temperature on the heated surface.';
       'Second input box is for the initial temperature of the object';
       'click on UPLOAD before using this App';
       'select the function use otherwise this App could not work';};
helpdlg(txt,dlgname);


% --------------------------------------------------------------------
function AppHelp_Callback(hObject, eventdata, handles)
% hObject    handle to AppHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function thermaldelta_Callback(hObject, eventdata, handles)
% hObject    handle to thermaldelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
input = str2num(get(hObject,'String'));
if(isempty(input))
    set(hObject,'String',0)
end
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of thermaldelta as text
%        str2double(get(hObject,'String')) returns contents of thermaldelta as a double


% --- Executes during object creation, after setting all properties.
function thermaldelta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thermaldelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function densitydelta_Callback(hObject, eventdata, handles)
% hObject    handle to densitydelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
input = str2num(get(hObject,'String'));
if(isempty(input))
    set(hObject,'String',0)
end
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of densitydelta as text
%        str2double(get(hObject,'String')) returns contents of densitydelta as a double


% --- Executes during object creation, after setting all properties.
function densitydelta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to densitydelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function heatdelta_Callback(hObject, eventdata, handles)
% hObject    handle to heatdelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
input = str2num(get(hObject,'String'));
if(isempty(input))
    set(hObject,'String',0)
end
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of heatdelta as text
%        str2double(get(hObject,'String')) returns contents of heatdelta as a double


% --- Executes during object creation, after setting all properties.
function heatdelta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to heatdelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function func5(handles)
    global a b c d e
%     x = get(handles.diyigeshu,'String');
%     y = get(handles.diergeshu,'String');
%     a = str2double(x);
%    
%     b = str2double(y);
%     c = a+b
%     d = num2str(c);
% %     set(handles.jieguo,'String',d);
%     guidata(hObject.handles);
%  a = str2double(get(handles.diyigeshu,'String'));
% b =str2double(get(handles.diergeshu,'String'));
A = createpde('thermal','steadystate');
% point1 = [3 4 -.075 .075 .075 -.075  -.1125 -.1125 .1125 .1125];
% point2 = [3 4 -.075 -.0075 -.0075 -.075  -.0975 -.0975 .0975 .0975];
% point3 = [3 4 .075 .0075 .0075 .075  -.0975 -.0975 .0975 .0975];

 point1 = [3 4 -.075 .075 .075 -.075  -.1125 -.1125 .1125 .1125];
 point2 = [3 4 -.075 -.0075 -.0075 -.075  -.0975 -.0975 .0975 .0975];
 point3 = [3 4 .075 .0075 .0075 .075  -.0975 -.0975 .0975 .0975];

reverse = [point1; point2;point3]';
g = decsg(reverse,'R1-R2-R3',['R1'; 'R2';'R3']');
geometryFromEdges(A,g);
 
 
thermalBC(A,'Edge',2,'HeatFlux',-22.2);
thermalBC(A,'Edge',3,'HeatFlux',-33.3);
thermalBC(A,'Edge',4,'HeatFlux',-33.3);
thermalBC(A,'Edge',5,'HeatFlux',-2.22);
thermalBC(A,'Edge',6,'HeatFlux',-2.22);
thermalBC(A,'Edge',7,'HeatFlux',-21.09);
thermalBC(A,'Edge',8,'HeatFlux',-21.09);
thermalBC(A,'Edge',9,'HeatFlux',-21.09);
thermalBC(A,'Edge',10,'HeatFlux',-21.09);
thermalBC(A,'Edge',11,'HeatFlux',-2.22);
thermalBC(A,'Edge',12,'HeatFlux',-2.22);
thermalBC(A,'Edge',1,'Temperature',a);
 
thermalProperties(A,'ThermalConductivity',e);
 
 
generateMesh(A,'Hmax',0.2);
 
R = solve(A);
T = R.Temperature;
figure
pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet');
axis equal
title 'Temperature, Steady State Solution'
 
 
B = createpde('thermal','transient');
 
point1 = [3 4 -.075 .075 .075 -.075  -.1125 -.1125 .1125 .1125];
point2 = [3 4 -.075 -.0075 -.0075 -.075  -.0975 -.0975 .0975 .0975];
point3 = [3 4 .075 .0075 .0075 .075  -.0975 -.0975 .0975 .0975];
gdm = [point1; point2;point3]';
g = decsg(gdm,'R1-R2-R3',['R1'; 'R2';'R3']');
geometryFromEdges(B,g);
 
 
thermalProperties(B,'MassDensity', c,'SpecificHeat',d, 'ThermalConductivity',e);
                             
                                
  
thermalBC(B,'Edge',2,'HeatFlux',-22.2);
thermalBC(B,'Edge',3,'HeatFlux',-33.3);
thermalBC(B,'Edge',4,'HeatFlux',-33.3);
thermalBC(B,'Edge',5,'HeatFlux',-2.22);
thermalBC(B,'Edge',6,'HeatFlux',-2.22);
thermalBC(B,'Edge',7,'HeatFlux',-21.09);
thermalBC(B,'Edge',8,'HeatFlux',-21.09);
thermalBC(B,'Edge',9,'HeatFlux',-21.09);
thermalBC(B,'Edge',10,'HeatFlux',-21.09);
thermalBC(B,'Edge',11,'HeatFlux',-2.22);
thermalBC(B,'Edge',12,'HeatFlux',-2.22);
thermalBC(B,'Edge',1,'Temperature',a);
 
msh = generateMesh(B,'Hmax',0.02);
 
tlist = 0:.1:50;
thermalIC(B,b);
R = solve(B,tlist);
T = R.Temperature;
 
 
getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
 
[~,nid] = getClosestNode( msh.Nodes, .5, 0 );
 
 
h = figure;
h.Position = [1 1 2 1].*h.Position;
subplot(1,2,1);
axis equal
pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet');
axis equal
title 'transient solution'
subplot(1,2,2);
axis equal
plot(tlist, T(nid,:));
grid on
title 'Temperature as different Time C';
xlabel 'seconds'
ylabel 'Temperature C'
close all
F(length(tlist)) = struct('cdata',[],'colormap',[]);
for tt=1:length(tlist)
pdeplot(A,'XYData',T(:,tt),'Contour','on','ColorMap','jet');
drawnow;
F(tt) = getframe;
end

function func6(handles)
global a b c d e

A = createpde('thermal','steadystate');
rect = [2 3 0 30 15 0 0 25.98];
tri1 = [2 3 1 9 5 1 1 7.93];
tri2 = [2 3 11 19 15 1 1 7.93];
tri3 = [2 3 21 29 25 1 1 7.93];
tri4 = [2 3 6 10 14 7.66 1 7.66];
tri5 = [2 3 16 20 24 7.66 1 7.66];
tri6 = [2 3 5 25 15 8.66 8.66 26];
gdm = [rect;tri1;tri2;tri3;tri4;tri5;tri6]';
g = decsg(gdm,'rect - tri1 - tri2 - tri3 - tri4 - tri5 - tri6',['rect';'tri1' ;'tri2';'tri3';'tri4';'tri5';'tri6']');
geometryFromEdges(A,g);


thermalBC(A,'Edge',2,'HeatFlux',-10);
thermalBC(A,'Edge',1,'Temperature',a);

thermalProperties(A,'ThermalConductivity',c);


generateMesh(A,'Hmax',0.2);

R = solve(A);
T = R.Temperature;
figure
pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet'); 
axis equal
title 'Temperature, Steady State Solution'
B = createpde('thermal','transient');

rect = [2 3 0 30 15 0 0 25.98];
tri1 = [2 3 1 9 5 1 1 7.93];
tri2 = [2 3 11 19 15 1 1 7.93];
tri3 = [2 3 21 29 25 1 1 7.93];
tri4 = [2 3 6 10 14 7.66 1 7.66];
tri5 = [2 3 16 20 24 7.66 1 7.66];
tri6 = [2 3 5 25 15 8.66 8.66 26];
gdm = [rect;tri1;tri2;tri3;tri4;tri5;tri6]';
g = decsg(gdm,'rect - tri1 - tri2 - tri3 - tri4 - tri5 - tri6',['rect';'tri1' ;'tri2';'tri3';'tri4';'tri5';'tri6']');
geometryFromEdges(B,g);

thermalProperties(B,'ThermalConductivity',c,'MassDensity',d,'SpecificHeat',e);

thermalBC(B,'Edge',1,'HeatFlux',-10);
% thermalBC(B,'Edge',6,'Temperature',@transientBCHeatedBlock);
thermalBC(B,'Edge',1,'Temperature',a);
msh = generateMesh(B,'Hmax',0.2);
tlist = 0:1000:400000;
thermalIC(B,b);
R = solve(B,tlist);
T = R.Temperature;


getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);

[~,nid] = getClosestNode( msh.Nodes, .5, 0 );


h = figure;
h.Position = [1 1 2 1].*h.Position;
subplot(1,2,1); 
axis equal  
pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet'); 
axis equal
title 'Temperature, Final Time, Transient Solution'
subplot(1,2,2); 
axis equal
plot(tlist, T(nid,:)); 
grid on
title 'Temperature at Right Edge as a Function of Time';
xlabel 'Time, seconds'
ylabel 'Temperature, degrees-Celsius'









%  
%  
%  
% k = @(~,state) 0.3+0.003*state.u;
%  
% thermalProperties(A,'ThermalConductivity',k);
%  
% R = solve(A);
% T = R.Temperature;
% figure
% pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet'); 
% axis equal
% title 'Temperature, Steady State Solution'
%  
%  
% thermalProperties(B,'ThermalConductivity',k,'MassDensity',1,'SpecificHeat',1);
%  
% thermalIC(B,0);
% R = solve(B,tlist);
% T = R.Temperature;
%  
%  
% h = figure;
% h.Position = [1 1 2 1].*h.Position;
% subplot(1,2,1);
% axis equal
% pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet');
% axis equal
% title 'Temperature, Final Time, Transient Solution'
% subplot(1,2,2);
% axis equal
% plot(tlist(1:size(T,2)), T(nid,:));
% grid on
% title 'Temperature at Right Edge as a Function of Time (Nonlinear)';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'
%  
%  
%  
% for animation
% pause(300)
% close all
% F(length(tlist)) = struct('cdata',[],'colormap',[]);
% for tt=1:length(tlist)
% pdeplot(A,'XYData',T(:,tt),'Contour','on','ColorMap','jet');
% drawnow;
% F(tt) = getframe;
% 
% end
% 
% 
% function func2(handles)
% A = createpde('thermal');
% g = decsg([3 4 -1.5 1.5 1.5 -1.5 0 0 .2 .2]');
% geometryFromEdges(A,g);
% 
% con = 51.5; 
% den = 7800; 
% capa = 482.5; 
% supp = 20000; 
% FCON = @(region,state) con*region.y;
% FCAPA = @(region,state) capa*region.y;
% FSUPP = @(region,state) supp*region.y;
% 
% thermalProperties(A,'ThermalConductivity',FCON);
% internalHeatSource(A,FSUPP);
% 
% 
% thermalBC(A,'Edge',2,'Temperature',100);
% 
% coff = @(region,~) 50*region.y;
% thermalBC(A,'Edge',3,'ConvectionCoefficient',coff,'AmbientTemperature',100);
%                    
% leftHF = @(region,~) 5000*region.y;
% thermalBC(A,'Edge',4,'HeatFlux',leftHF);
% 
% 
% generateMesh(A,'Hmax',0.2);
% 
% result = solve(A);
% T = result.Temperature;
% figure; 
% pdeplot(A,'XYData',T,'Contour','on'); 
% axis equal
% title 'Steady State Temperature';
% 
% B = createpde('thermal','transient');
% 
% g = decsg([3 4 -1.5 1.5 1.5 -1.5 0 0 .2 .2]');
% geometryFromEdges(B,g);
% 
% generateMesh(B,'Hmax',0.2);
% 
% thermalProperties('SpecificHeat',FCAPA ,B,'ThermalConductivity',FCON,'MassDensity',den);                           
% 
% internalHeatSource(B,FSUPP);
% 
% thermalBC(B,'Edge',2,'Temperature',100);
% thermalBC(B,'Edge',3,'ConvectionCoefficient',coff,'AmbientTemperature',100);
% thermalBC(B,'Edge',4,'HeatFlux',leftHF);
% 
% tfinal = 20000;
% tlist = 0:100:tfinal;
% thermalIC(B,0);
% B.SolverOptions.ReportStatistics = 'on';
% 
% result = solve(B,tlist);
% T = result.Temperature;
% 
% figure; 
% pdeplot(B,'XYData',T(:,end),'Contour','on'); 
% axis equal
% title(sprintf('Transient Temperature at Final Time (%g seconds)',tfinal));
% 
% p = B.Mesh.Nodes;
% nodesLeftEnd  = find(p(1,:) < -1.5+eps);
% nodeCenter = nodesLeftEnd(p(2,nodesLeftEnd) < eps);
% nodeOuter = nodesLeftEnd(p(2,nodesLeftEnd) > 0.2-eps);
% 
% figure; 
% plot(tlist,T(nodeCenter,:)); 
% hold all 
% plot(tlist,T(nodeOuter,:),'--');
% title 'Temperature at Left End as a Function of Time'
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-C'
% grid on;
% legend('Center Axis','Outer Surface','Location','SouthEast');
% 
% 
% function func3(handles);
% 
% A = createpde('thermal','steadystate');
% g = decsg([3 4 -0.5 0.5 0.5 -0.5 -4.5 -4.5 4.5 4.5]');
% geometryFromEdges(A,g);
%  
% 
% 
%  
% thermalBC(A,'Edge',2,'HeatFlux',-22.2);
% thermalBC(A,'Edge',3,'HeatFlux',-33.3);
% thermalBC(A,'Edge',4,'HeatFlux',-33.3);
% thermalBC(A,'Edge',5,'HeatFlux',-2.22);
% thermalBC(A,'Edge',6,'HeatFlux',-2.22);
% thermalBC(A,'Edge',7,'HeatFlux',-21.09);
% thermalBC(A,'Edge',8,'HeatFlux',-21.09);
% thermalBC(A,'Edge',9,'HeatFlux',-21.09);
% thermalBC(A,'Edge',10,'HeatFlux',-21.09);
% thermalBC(A,'Edge',11,'HeatFlux',-2.22);
% thermalBC(A,'Edge',12,'HeatFlux',-2.22);
% thermalBC(A,'Edge',1,'Temperature',400);
%  
% thermalProperties(A,'ThermalConductivity',51.5);
%  
%  
% generateMesh(A,'Hmax',0.2);
%  
% R = solve(A);
% T = R.Temperature;
% figure
% pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet');
% axis equal
% title 'Temperature, Steady State Solution'
%  
%  
% B = createpde('thermal','transient');
% 
% g = decsg([3 4 -0.5 0.5 0.5 -0.5 -4.5 -4.5 4.5 4.5]');
% 
% geometryFromEdges(B,g);
%  
%  
% thermalProperties(B,'MassDensity',7800, 'ThermalConductivity',51.5,'SpecificHeat',482.5);
%                              
%                                 
%  
%  
% thermalBC(B,'Edge',2,'HeatFlux',-10);
% thermalBC(B,'Edge',1,'Temperature',400);
%  
% msh = generateMesh(B,'Hmax',0.2);
%  
% tlist = 0:.1:10;
% thermalIC(B,0);
% R = solve(B,tlist);
% T = R.Temperature;
%  
%  
% getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
%  
% [~,nid] = getClosestNode( msh.Nodes, .5, 0 );
%  
%  
% h = figure;
% h.Position = [1 1 2 1].*h.Position;
% subplot(1,2,1);
% axis equal
% pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet');
% axis equal
% title 'Temperature, Final Time, Transient Solution'
% subplot(1,2,2);
% axis equal
% plot(tlist, T(nid,:));
% grid on
% title 'Temperature at Right Edge as a Function of Time';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'
%  
%  
%  
% k = @(~,state) 0.3+0.003*state.u;
%  
% thermalProperties(A,'ThermalConductivity',k);
%  
% R = solve(A);
% T = R.Temperature;
% figure
% pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet'); 
% axis equal
% title 'Temperature, Steady State Solution'
%  
%  
% thermalProperties(B,'ThermalConductivity',k,'MassDensity',1,'SpecificHeat',1);
%  
% thermalIC(B,b);
% R = solve(B,tlist);
% T = R.Temperature;
%  
%  
% h = figure;
% h.Position = [1 1 2 1].*h.Position;
% subplot(1,2,1);
% axis equal
% pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet');
% axis equal
% title 'Temperature, Final Time, Transient Solution'
% subplot(1,2,2);
% axis equal
% plot(tlist(1:size(T,2)), T(nid,:));
% grid on
% title 'Temperature at Right Edge as a Function of Time (Nonlinear)';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'
% F(length(tlist)) = struct('cdata',[],'colormap',[]);
% 
% 
% 
%  
%  
% 
% function func4(handles)
% 
% A = createpde('thermal','steadystate');
% point2 = [2 3 -8 0 8 1 9 1];
% point1 = [2 3 -10 0 10 0 10 0];
% 
% gdm = [point1; point2]';
% g = decsg(gdm,'R1-R2',['R1'; 'R2']');
% geometryFromEdges(A,g);
% 
% figure;
% pdegplot(A,'EdgeLabels','on'); 
% axis([-20 20 -20 20]);
% title 'Block Geometry With Edge Labels Displayed'
% 
% thermalBC(A,'Edge',2,'HeatFlux',-10);
% thermalBC(A,'Edge',1,'Temperature',400);
% 
% thermalProperties(A,'ThermalConductivity',1);
% 
% 
% generateMesh(A,'Hmax',0.2);
% 
% R = solve(A);
% T = R.Temperature;
% figure
% pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet'); 
% axis equal
% title 'Temperature, Steady State Solution'
% 
% 
% B = createpde('thermal','transient');
% 
% r1 = [2 3 -1 0 1 0 1 0];
% point2 = [2 3 -10 0 10 0 10 0];
% 
% gdm = [r1; point2]';
% g = decsg(gdm,'R2-R1',['R1'; 'R2']');
% geometryFromEdges(B,g);
% 
% 
% thermalProperties(B,'MassDensity',1, 'ThermalConductivity',1,'SpecificHeat',1);
%                               
%                                 
% 
% 
% thermalBC(B,'Edge',2,'HeatFlux',-10);
% thermalBC(B,'Edge',1,'Temperature',@transientBCHeatedBlock);
% 
% msh = generateMesh(B,'Hmax',0.2);
% 
% tlist = 0:.1:5;
% thermalIC(B,0);
% R = solve(B,tlist);
% T = R.Temperature;
% 
% 
% getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
% 
% [~,nid] = getClosestNode( msh.Nodes, .5, 0 );
% 
% 
% h = figure;
% h.Position = [1 1 2 1].*h.Position;
% subplot(1,2,1); 
% axis equal
% pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet'); 
% axis equal
% title 'Temperature, Final Time, Transient Solution'
% subplot(1,2,2); 
% axis equal
% plot(tlist, T(nid,:)); 
% grid on
% title 'Temperature at Right Edge as a Function of Time';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'
% 
% 
% 
% k = @(~,state) 0.3+0.003*state.u;
% 
% thermalProperties(A,'ThermalConductivity',k);
%  
% 
% R = solve(A);
% T = R.Temperature;
% figure
% pdeplot(A,'XYData',T,'Contour','on','ColorMap','jet');  
% axis equal
% title 'Temperature, Steady State Solution'
% 
% 
% thermalProperties(B,'ThermalConductivity',k,'MassDensity',1,'SpecificHeat',1);
%  
% 
% thermalIC(B,0);
% R = solve(B,tlist);
% T = R.Temperature;
% 
% 
% h = figure;
% h.Position = [1 1 2 1].*h.Position;
% subplot(1,2,1); 
% axis equal
% pdeplot(B,'XYData',T(:,end),'Contour','on','ColorMap','jet'); 
% axis equal
% title 'Temperature, Final Time, Transient Solution'
% subplot(1,2,2); 
% axis equal
% plot(tlist(1:size(T,2)), T(nid,:)); 
% grid on
% title 'Temperature at Right Edge as a Function of Time (Nonlinear)';
% xlabel 'Time, seconds'
% ylabel 'Temperature, degrees-Celsius'
% 
% --- Executes just before GUI1 is made visible.
% function GUI1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI1 (see VARARGIN)
% 
% Choose default command line output for GUI1
% handles.output = hObject;
% 
% Update handles structure
% guidata(hObject, handles);
% 
% UIWAIT makes GUI1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% 
% --- Outputs from this function are returned to the command line.
% function varargout = GUI1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% Get default command line output from handles structure
% varargout{1} = handles.output;
% 
% 