function varargout = FD_sim_gui(varargin)
% FD_SIM_GUI MATLAB code for FD_sim_gui.fig
%      FD_SIM_GUI, by itself, creates a new FD_SIM_GUI or raises the existing
%      singleton*.
%
%      H = FD_SIM_GUI returns the handle to a new FD_SIM_GUI or the handle to
%      the existing singleton*.
%
%      FD_SIM_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FD_SIM_GUI.M with the given input arguments.
%
%      FD_SIM_GUI('Property','Value',...) creates a new FD_SIM_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FD_sim_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FD_sim_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FD_sim_gui

% Last Modified by GUIDE v2.5 05-May-2020 13:50:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FD_sim_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @FD_sim_gui_OutputFcn, ...
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


% --- Executes just before FD_sim_gui is made visible.
function FD_sim_gui_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FD_sim_gui (see VARARGIN)

% Obtain input from FD_seis_sim
handles.output = hObject;
handles.Vp=(varargin{1});
handles.Vs=(varargin{2});
handles.rho=(varargin{3});
handles.X=(varargin{4});
handles.Y=(varargin{5});
handles.Z=(varargin{6});
handles.topoload=(varargin{7});
handles.att_load=(varargin{8});

% Obtain topography (if given as input), or check if dummy values are given
if handles.topoload==1 && isequal(length(varargin{9}),length(varargin{10}),length(varargin{11}),1)==0
    handles.topo=(varargin{9});
    handles.X_topo=(varargin{10});
    handles.Y_topo=(varargin{11});
else
     handles.topoload=0;
     handles.topo=0;
     handles.X_topo=0;
     handles.Y_topo=0;
end


% Check if attenuation values are given as input
handles.att_type=1; % Default case
if handles.att_load==1  %(case 1: only Q0 or TP is given)
    handles.Q_tsl=varargin{12};
    set(handles.edit_Q_tsl,'enable','off')
    handles.att_type=1;
else % set some default value
    handles.Q_tsl=10;
end

if handles.att_load==2  %(case 2: Q0/f0 or TP/TS are given)
    handles.Q_tsl=varargin{12};
    set(handles.edit_Q_tsl,'enable','off')
    handles.f_tpp=varargin{13};
    set(handles.edit_f_tpp,'enable','off')
    handles.att_type=1;
else % set some default value
    handles.f_tpp=0.5;
end

if handles.att_load==3 %(case 3: TP/TS and tsl are given)
    handles.Q_tsl=varargin{12};
    set(handles.edit_Q_tsl,'enable','off')
    handles.f_tpp=varargin{13};
    set(handles.edit_f_tpp,'enable','off')
    handles.Q_tsl2=varargin{14};    
    handles.att_type=2;
    set(handles.edit_q_tsl2,'enable','off')
    set(handles.Button_att1,'enable','off')
    set(handles.Button_att2,'enable','off')
else % set some default value
    handles.Q_tsl2=0.22;
end

% Default parameters
handles.source=[90 45 270 1 18];    % Strike, dip, rake, scalar moment, scale
handles.fault=[20 10 2.5]; % fault length, width, rupture velocity
handles.dur=0.1; % initial pulse sigma (for gausisan) or period (for Ricker)
handles.maxt=2000; % Number of time-steps to simulate
handles.dt=0.01; % sampling period
handles.posx=round(length(handles.Vp(1,:,1))./2); % Set middle of simulation space as source location
handles.posy=round(length(handles.Vp(:,1,1))./2); % Set middle of simulation space as source location
handles.posz=round(length(handles.Vp(1,1,:))./2); % Set middle of simulation space as source location
handles.xrup=0.5;handles.zrup=0.5; % initial X and Y (Z) location of rupture initiation within fault (0.5=middle of fault) 
handles.gpumode=1; % Gpu processing on
handles.scheme=1; % Spatial partial derivative scheme (1= 2nd order, 2=4th order)
handles.stype=1; % Source type (1=point, 2=fault)
handles.outpath=pwd; % initial output path
handles.savemode=2; % 1=save whole 3D wavefield, 2=save only surface values
handles.savetype=1; % 1=save as .mat, 2=save as binary
handles.skipt=10; % save date every 10 simulated timesteps
handles.N_abc=10; % Number of nodes to use as transition zone in the side boundaries 
handles.b_abc=0.5; % 
guidata(hObject,handles)


% Set up an initial default source pulse
fs=1000; %sampling frequency
t=-1:1/fs:1; %time base 
variance=handles.dur^2; 
x=1/(sqrt(2*pi*variance))*(exp(-t.^2/(2*variance)));
x(abs(x)<0.000001)=0;
temp=find(x,1);
handles.pulse(:,1)=t+abs(t(temp(1)));
handles.pulse(:,2)=x./max(abs(x));
handles.pulse(:,2)=x./max(abs(x));
handles.stf_type=1;

% Insert initial values into GUI
set(handles.edit_outpath,'String',pwd)  
set(handles.edit_mt1,'String',num2str(handles.source(1)))
set(handles.edit_mt2,'String',num2str(handles.source(2)))
set(handles.edit_mt3,'String',num2str(handles.source(3)))
set(handles.edit_mt4,'String',num2str(handles.source(4)))
set(handles.edit_mt5,'String',num2str(handles.source(5)))

set(handles.edit_fault1,'String',num2str(handles.fault(1)))
set(handles.edit_fault2,'String',num2str(handles.fault(2)))
set(handles.edit_fault3,'String',num2str(handles.fault(3)))

set(handles.edit_stf,'String',num2str(handles.dur))
set(handles.edit_maxt,'String',num2str(handles.maxt))
set(handles.edit_dt,'String',num2str(handles.dt))

set(handles.edit_ix,'String',num2str(handles.posx))
set(handles.edit_iy,'String',num2str(handles.posy))
set(handles.edit_iz,'String',num2str(handles.posz))

set(handles.edit_N_abc,'String',num2str(handles.N_abc))
set(handles.edit_b_abc,'String',num2str(handles.b_abc))

set(handles.edit_txx,'enable','off')
set(handles.edit_tyy,'enable','off')
set(handles.edit_tzz,'enable','off')
set(handles.edit_txy,'enable','off')
set(handles.edit_txz,'enable','off')
set(handles.edit_tyz,'enable','off')
set(handles.edit_fault1,'enable','off')
set(handles.edit_fault2,'enable','off')
set(handles.edit_fault3,'enable','off')
set(handles.edit_xrup,'enable','off')
set(handles.edit_zrup,'enable','off')
set(handles.edit_q_tsl2,'enable','off')

% Update handles structure
guidata(hObject, handles);




% --- Outputs from this function are returned to the command line.
function varargout = FD_sim_gui_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when selected object is changed in Group1_stype.
function Group1_stype_SelectionChangedFcn(hObject, ~, handles)
% hObject    handle to the selected object in Group1_stype 
% handles    structure with handles and user data (see GUIDATA)
if (hObject == handles.Button_stype1)
    handles.stype=1;
    set(handles.edit_fault1,'enable','off')
    set(handles.edit_fault2,'enable','off')
    set(handles.edit_fault3,'enable','off')
    set(handles.edit_xrup,'enable','off')
    set(handles.edit_zrup,'enable','off')
elseif (hObject == handles.Button_stype2)
    handles.stype=2;
    set(handles.edit_fault1,'enable','on')
    set(handles.edit_fault2,'enable','on')
    set(handles.edit_fault3,'enable','on')
    set(handles.edit_xrup,'enable','on')
    set(handles.edit_zrup,'enable','on')
end
guidata(hObject,handles)
 
   

function edit_mt1_Callback(hObject, ~, handles)
% hObject    handle to edit_mt1 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit Strike
val=get(handles.edit_mt1,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.source(1)));
    errordlg('Input must be a number between -360 and 360','Error');
    return
elseif abs(str2double(val))>360
    set(hObject, 'String',num2str(handles.source(1)));
    errordlg('Input must be a number between -360 and 360','Error');
    return
end
handles.source(1)=str2double(val);
[txx,txy,txz,tyy,tyz,tzz]=sdr2mt(handles.source(1),handles.source(2),handles.source(3),handles.source(4));
set(handles.edit_txx,'String',num2str(txx))
set(handles.edit_txy,'String',num2str(txy))
set(handles.edit_txz,'String',num2str(txz))
set(handles.edit_tyy,'String',num2str(tyy))
set(handles.edit_tyz,'String',num2str(tyz))
set(handles.edit_tzz,'String',num2str(tzz))

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_mt1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_mt1 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mt2_Callback(hObject, ~, handles)
% hObject    handle to edit_mt2 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit dip
val=get(handles.edit_mt2,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.source(2)));
    errordlg('Input must be a number between -90 and 90','Error');
    return
elseif abs(str2double(val))>90
    set(hObject, 'String',num2str(handles.source(2)));
    errordlg('Input must be a number between -90 and 90','Error');
    return
end
handles.source(2)=str2double(val);
[txx,txy,txz,tyy,tyz,tzz]=sdr2mt(handles.source(1),handles.source(2),handles.source(3),handles.source(4));
set(handles.edit_txx,'String',num2str(txx))
set(handles.edit_txy,'String',num2str(txy))
set(handles.edit_txz,'String',num2str(txz))
set(handles.edit_tyy,'String',num2str(tyy))
set(handles.edit_tyz,'String',num2str(tyz))
set(handles.edit_tzz,'String',num2str(tzz))

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_mt2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_mt2 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mt3_Callback(hObject, ~, handles)
% hObject    handle to edit_mt3 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit rake
val=get(handles.edit_mt3,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.source(3)));
    errordlg('Input must be a number between -360 and 360','Error');
    return
elseif abs(str2double(val))>360
    set(hObject, 'String',num2str(handles.source(3)));
    errordlg('Input must be a number between -360 and 360','Error');
    return
end
handles.source(3)=str2double(val);
[txx,txy,txz,tyy,tyz,tzz]=sdr2mt(handles.source(1),handles.source(2),handles.source(3),handles.source(4));
set(handles.edit_txx,'String',num2str(txx))
set(handles.edit_txy,'String',num2str(txy))
set(handles.edit_txz,'String',num2str(txz))
set(handles.edit_tyy,'String',num2str(tyy))
set(handles.edit_tyz,'String',num2str(tyz))
set(handles.edit_tzz,'String',num2str(tzz))
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function edit_mt3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_mt3 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mt4_Callback(hObject, ~, handles)
% hObject    handle to edit_mt4 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit scalar moment (without exponent)
val=get(handles.edit_mt4,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.source(4)));
    errordlg('Input must be a number','Error');
end
handles.source(4)=str2double(val);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_mt4_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_mt4 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mt5_Callback(hObject, ~, handles)
% hObject    handle to edit_mt5 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit moment scale (exponent)
val=get(handles.edit_mt5,'String');
if isnan(val)
    set(hObject, 'String',num2str(handles.source(5)));
    errordlg('Input must be a number','Error');
end
handles.source(5)=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_mt5_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_mt5 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_txx_Callback(~, ~, ~)
% hObject    handle to edit_txx (see GCBO)
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function edit_txx_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_txx (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tyy_Callback(~, ~, ~)
% hObject    handle to edit_tyy (see GCBO)
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function edit_tyy_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_tyy (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tzz_Callback(~, ~, ~)
% hObject    handle to edit_tzz (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit_tzz_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_tzz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_txy_Callback(~, ~, ~)
% hObject    handle to edit_txy (see GCBO)
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit_txy_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_txy (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_txz_Callback(~, ~, ~)
% hObject    handle to edit_txz (see GCBO)
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit_txz_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_txz (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tyz_Callback(~, ~, ~)
% hObject    handle to edit_tyz (see GCBO)
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit_tyz_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_tyz (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fault1_Callback(hObject, ~, handles)
% hObject    handle to edit_fault1 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit fault length
val=get(handles.edit_fault1,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.fault(1)));
    errordlg('Input must be a number','Error');
end
handles.fault(1)=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_fault1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_fault1 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fault2_Callback(hObject, ~, handles)
% hObject    handle to edit_fault2 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit fault width
val=get(handles.edit_fault2,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.fault(2)));
    errordlg('Input must be a number','Error');
end
handles.fault(2)=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_fault2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_fault2 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fault3_Callback(hObject, ~, handles)
% hObject    handle to edit_fault3 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit fault minimum depth from surface
val=get(handles.edit_fault3,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.fault(3)));
    errordlg('Input must be a number','Error');
end
handles.fault(3)=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_fault3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_fault3 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, ~, handles)
% hObject    handle to popupmenu1 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Select input source time function type and plot based on input parameters
% (sigma or period)
handles.stf_type=get(hObject,'Value');
switch get(hObject,'Value') 
    case 1 % Gaussian
        fs=1000; %sampling frequency
        sigma=handles.dur;
        t=-1:1/fs:1; %time base 
        variance=sigma^2;
        x=1/(sqrt(2*pi*variance))*(exp(-t.^2/(2*variance)));
        [x2, index] = unique(x(1:1001));
        a=interp1(x2./max(x2),t(index),0.5);

        figure
        fill([0 -a+1 -a+1 a+1 a+1 0],[0 0 0.5 0.5 0 0],'r','FaceAlpha',0.2)
        hold on
        plot(t+1,x./max(x),'linewidth',1)
        plot([a+1 -a+1], [0.5 0.5],'x')
        plot([a+1 -a+1], [0.5 0.5],'r')
        text(1,0.7,['Width: ' num2str(abs(2*a)) ' sec'],'HorizontalAlignment','center')
        xlabel('Time (s)')
        ylabel('Normalized Amplitude')
        title('Gaussian source pulse')
    case 2 % Gaussian 1st Derivative
        fs=1000; %sampling frequency
        sigma=handles.dur;
        t=-1:1/fs:1; %time base 
        variance=sigma^2;
        x=1/(sqrt(2*pi*variance))*(exp(-t.^2/(2*variance)));
        [x2, index] = unique(x(1:1001));
        a=interp1(x2./max(x2),t(index),0.5);
        x=diff(x);
        x(end+1)=x(end);
        
        figure
        fill([0 -a+1 -a+1 a+1 a+1 0],[-1 -1 1 1 -1 -1],'r','FaceAlpha',0.2)
        hold on
        plot(t+1,x./max(x),'linewidth',1)
        plot([a+1 -a+1], [0.5 0.5],'x')
        plot([a+1 -a+1], [0.5 0.5],'r')
        text(1,0.7,['Width: ' num2str(abs(2*a)) ' sec'],'HorizontalAlignment','center')
        xlabel('Time (s)')
        ylabel('Normalized Amplitude')
        title('Gaussian derivative source pulse')
    case 3 % Ricker
        fs=1000;
        tp=handles.dur;
        ts=0;
        t=-1:1/fs:1; %time base
        a=(pi.*((t-ts)./tp)).^2;
        x=(sqrt(pi)/2).*(a-0.5).*(exp(-a));
        figure
        hold on
        plot(t+1,x./max(abs(x)),'linewidth',1)
        text(0.5,-0.7,['Dominant frequency: ' num2str(1/tp) ' Hz'],'HorizontalAlignment','center')
        xlabel('Time (s)')
        ylabel('Normalized Amplitude')
        title('Ricker source pulse')
end

x(abs(x)<0.000001)=0;
temp=find(x,1);
handles.pulse(:,1)=t+abs(t(temp(1)));
handles.pulse(:,2)=x./max(abs(x));
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stf_Callback(hObject, ~, handles)
% hObject    handle to edit_stf (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit sigma for Gaussian source time function of Period for Ricker
val=get(handles.edit_stf,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.dur));
    errordlg('Input must be a number','Error');
    return
end
handles.dur=str2double(val);

if handles.stf_type==1
    fs=1000; %sampling frequency
    sigma=handles.dur;
    t=-1:1/fs:1; %time base 
    variance=sigma^2;
    x=1/(sqrt(2*pi*variance))*(exp(-t.^2/(2*variance)));
    [x2, index] = unique(x(1:1001));
    a=interp1(x2./max(x2),t(index),0.5);

    figure
    fill([0 -a+1 -a+1 a+1 a+1 0],[0 0 0.5 0.5 0 0],'r','FaceAlpha',0.2)
    hold on
    plot(t+1,x./max(x),'linewidth',1)
    plot([a+1 -a+1], [0.5 0.5],'x')
    plot([a+1 -a+1], [0.5 0.5],'r')
    text(1,0.7,['Width: ' num2str(abs(2*a)) ' sec'],'HorizontalAlignment','center')
    xlabel('Time (s)')
    ylabel('Normalized Amplitude')
    title('Gaussian source pulse')
elseif handles.stf_type==2
    fs=1000; %sampling frequency
    sigma=handles.dur;
    t=-1:1/fs:1; %time base 
    variance=sigma^2;
    x=1/(sqrt(2*pi*variance))*(exp(-t.^2/(2*variance)));
    [x2, index] = unique(x(1:1001));
    a=interp1(x2./max(x2),t(index),0.5);
    x=diff(x);
    x(end+1)=x(end);
        
    figure
    fill([0 -a+1 -a+1 a+1 a+1 0],[-1 -1 1 1 -1 -1],'r','FaceAlpha',0.2)
    hold on
    plot(t+1,x./max(x),'linewidth',1)
    plot([a+1 -a+1], [0.5 0.5],'x')
    plot([a+1 -a+1], [0.5 0.5],'r')
    text(1,0.7,['Width: ' num2str(abs(2*a)) ' sec'],'HorizontalAlignment','center')
    xlabel('Time (s)')
    ylabel('Normalized Amplitude')
    title('Gaussian derivative source pulse')
elseif handles.stf_type==3
    fs=1000;
    tp=handles.dur;
    ts=0;
    t=-1:1/fs:1; %time base
    a=(pi.*((t-ts)./tp)).^2;
    x=(sqrt(pi)/2).*(a-0.5).*(exp(-a));
    figure
    hold on
    plot(t+1,x./max(abs(x)),'linewidth',1)
    text(0.5,-0.7,['Dominant frequency: ' num2str(1/tp) ' Hz'],'HorizontalAlignment','center')
    xlabel('Time (s)')
    ylabel('Normalized Amplitude')
    title('Ricker source pulse')
end

x(abs(x)<0.000001)=0;
temp=find(x,1);
handles.pulse(:,1)=t+abs(t(temp(1)));
handles.pulse(:,2)=x./max(abs(x));
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function edit_stf_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_stf (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxt_Callback(hObject, ~, handles)
% hObject    handle to edit_maxt (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit number of simulation steps
val=get(handles.edit_maxt,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.maxt));
    errordlg('Input must be a number','Error');
    return
end
handles.maxt=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_maxt_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_maxt (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dt_Callback(hObject, ~, handles)
% hObject    handle to edit_dt (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit sampling period
val=get(handles.edit_dt,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.dt));
    errordlg('Input must be a number','Error');
    return
end
handles.dt=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_dt_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_dt (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ix_Callback(hObject, ~, handles)
% hObject    handle to edit_ix (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit source X node index
val=get(handles.edit_ix,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.posx));
    errordlg('Input must be a number','Error');
    return
elseif str2num(val)>length(handles.Vp(1,:,1)) || str2num(val)<1
    set(hObject, 'String',num2str(handles.posx));
    errordlg('Source node index must be 0 < index < Nx','Error');
    return
end
handles.posx=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_ix_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_ix (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_iy_Callback(hObject, ~, handles)
% hObject    handle to edit_iy (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit source Y node index
val=get(handles.edit_iy,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.posy));
    errordlg('Input must be a number','Error');
    return
elseif str2num(val)>length(handles.Vp(:,1,1)) || str2num(val)<1
    set(hObject, 'String',num2str(handles.posy));
    errordlg('Source node index must be 0 < index < Ny','Error');
    return
end
handles.posy=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_iy_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_iy (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_iz_Callback(hObject, ~, handles)
% hObject    handle to edit_iz (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit source Z node index
val=get(handles.edit_iz,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.posz));
    errordlg('Input must be a number','Error');
    return
elseif str2num(val)>length(handles.Vp(1,1,:)) || str2num(val)<1
    set(hObject, 'String',num2str(handles.posz));
    errordlg('Source node index must be 0 < index < Ny','Error');
    return
end
handles.posz=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_iz_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_iz (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xrup_Callback(hObject, ~, handles)
% hObject    handle to edit_xrup (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit rupture X initiation point within fault
val=get(handles.edit_xrup,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.xrup));
    errordlg('Input must be a number','Error');
    return
elseif str2num(val)>1 || str2num(val)<0
    set(hObject, 'String',num2str(handles.xrup));
    errordlg('rupture source must be 0 < s_rupt <1 ','Error');
    return
end
handles.xrup=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_xrup_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_xrup (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_zrup_Callback(hObject, ~, handles)
% hObject    handle to edit_zrup (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit rupture Y (Z) initiation point within fault
val=get(handles.edit_zrup,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.zrup));
    errordlg('Input must be a number','Error');
    return
elseif str2num(val)>1 || str2num(val)<0
    set(hObject, 'String',num2str(handles.zrup));
    errordlg('rupture source must be 0 < s_rupt <1 ','Error');
    return
end
handles.zrup=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_zrup_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_zrup (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup3.
function uibuttongroup3_SelectionChangedFcn(hObject, ~, handles)
% hObject    handle to the selected object in uibuttongroup3 
% handles    structure with handles and user data (see GUIDATA)
% Select attenuation type
if (hObject == handles.Button_att1)
    handles.att_type=1;
    set(handles.text_att1,'String','Q0')
    set(handles.text_att2,'String','f0')
    set(handles.edit_q_tsl2,'enable','off')
elseif (hObject == handles.Button_att2)
    handles.att_type=2;
    set(handles.text_att1,'String','TP')
    set(handles.text_att2,'String','TS')
    set(handles.edit_q_tsl2,'enable','on')
end
guidata(hObject,handles)



function edit_Q_tsl_Callback(hObject, ~, handles)
% hObject    handle to edit_Q_tsl (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit Q0 (if attenuation type =1) or TP (if attenuation type =2)
val=get(handles.edit_Q_tsl,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.Q_tsl));
    errordlg('Input must be a number','Error');
    return
end
handles.Q_tsl=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_Q_tsl_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_Q_tsl (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_f_tpp_Callback(hObject, ~, handles)
% hObject    handle to edit_f_tpp (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Edit f0 (if attenuation type =1) or TS (if attenuation type =2)
val=get(handles.edit_f_tpp,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.f_tpp));
    errordlg('Input must be a number','Error');
    return
end
handles.f_tpp=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_f_tpp_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_f_tpp (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_q_tsl2_Callback(hObject, ~, handles)
% hObject    handle to edit_q_tsl2 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
val=get(handles.edit_q_tsl2,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.Q_tsl2));
    errordlg('Input must be a number','Error');
    return
end
handles.Q_tsl2=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_q_tsl2_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_q_tsl2 (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_gpu.
function checkbox_gpu_Callback(hObject, ~, handles)
% hObject    handle to checkbox_gpu (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Select/Deselect GPU processing
val=get(handles.checkbox_gpu,'Value');
if val==1
    handles.gpumode=1;
else
    handles.gpumode=0;
end
guidata(hObject,handles)


% --- Executes when selected object is changed in Group1_scheme.
function Group1_scheme_SelectionChangedFcn(hObject, ~, handles)
% hObject    handle to the selected object in Group1_scheme 
% handles    structure with handles and user data (see GUIDATA)
% Select spatial derivative scheme
if (hObject == handles.Button_scheme1)
    handles.scheme=1;
elseif (hObject == handles.Button_scheme2)
    handles.scheme=2;
end
guidata(hObject,handles)



function edit_outpath_Callback(hObject, ~, handles)
% hObject    handle to edit_outpath (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Set output path
val=get(handles.edit_outpath,'String');
handles.outpath=val;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_outpath_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_outpath (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_path.
function button_path_Callback(hObject, ~, handles)
% hObject    handle to button_path (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Get output path
handles.outpath=uigetdir;
set(handles.edit_outpath,'String',handles.outpath)
guidata(hObject,handles)


% --- Executes on button press in checkbox_savemod.
function checkbox_savemod_Callback(hObject, ~, handles)
% hObject    handle to checkbox_savemod (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Save whole 3D or only surface Vx, Vy and Vz simulated values
val=get(handles.checkbox_savemod,'Value');
if val==1
    handles.savemode=2;
else
    handles.savemode=1;
end
guidata(hObject,handles)


% --- Executes on button press in checkbox_savetype.
function checkbox_savetype_Callback(hObject, ~, handles)
% hObject    handle to checkbox_savetype (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Set saving mode (mat or binary)
val=get(handles.checkbox_savetype,'Value');
if val==1
    handles.savetype=2;
else
    handles.savetype=1;
end
guidata(hObject,handles)



function edit_skipt_Callback(hObject, ~, handles)
% hObject    handle to edit_skipt (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Saving period
val=get(handles.edit_skipt,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.skipt));
    errordlg('Input must be a number','Error');
    return
end
handles.skipt=str2double(val);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_skipt_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_skipt (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_N_abc_Callback(hObject, ~, handles)
% hObject    handle to edit_N_abc (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Set the number of nodes to use as transition zone in the side boundaries
% in order to apply the absorbing boundary conditions
val=get(handles.edit_N_abc,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.N_abc));
    errordlg('Input must be a number','Error');
    return
end
handles.N_abc=str2double(val);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_N_abc_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_N_abc (see GCBO)
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_b_abc_Callback(hObject, ~, handles)
% hObject    handle to edit_b_abc (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
val=get(handles.edit_b_abc,'String');
if isempty(str2num(val))
    set(hObject, 'String',num2str(handles.b_abc));
    errordlg('Input must be a number','Error');
    return
end
handles.b_abc=str2double(val);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_b_abc_CreateFcn(hObject, ~, ~)
% hObject    handle to edit_b_abc (see GCBO)
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Button_close.
function Button_close_Callback(~, ~, ~)
% hObject    handle to Button_close (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Closes the GUI
close

% --- Executes on button press in Button_run.
function Button_run_Callback(~, ~, handles)
% hObject    handle to Button_run (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
% Calls the appropriate script on order to do the simulation based on the
% desired parameters

if handles.stype==1 && handles.att_type==1 && handles.scheme==1 
    % Point source, 2nd order spatial derivatives, constant Q
    FD_sim_point_2(handles.Vp,handles.Vs,handles.rho,...
        handles.X,handles.Y,handles.Z,handles.source,handles.pulse,...
        handles.Q_tsl,handles.f_tpp,handles.dt,handles.maxt,handles.posx,...
        handles.posy,handles.posz,handles.gpumode,handles.N_abc,handles.b_abc,...
        handles.topoload,handles.topo,handles.X_topo,handles.Y_topo,...
        handles.savemode,handles.savetype,handles.skipt,handles.outpath);
elseif handles.stype==1 && handles.att_type==2 && handles.scheme==1
    % Point source, 2nd order spatial derivatives, Q based on relaxation functions
    FD_sim_point_2_Q2(handles.Vp,handles.Vs,handles.rho,...
        handles.X,handles.Y,handles.Z,handles.source,handles.pulse,...
        handles.Q_tsl,handles.f_tpp,handles.Q_tsl2,handles.dt,handles.maxt,handles.posx,...
        handles.posy,handles.posz,handles.gpumode,handles.N_abc,handles.b_abc,...
        handles.topoload,handles.topo,handles.X_topo,handles.Y_topo,...
        handles.savemode,handles.savetype,handles.skipt,handles.outpath);
elseif handles.stype==1 && handles.att_type==2 && handles.scheme==2
    % Point source, 4th order spatial derivatives, Q based on relaxation functions
    FD_sim_point_4_Q2(handles.Vp,handles.Vs,handles.rho,...
        handles.X,handles.Y,handles.Z,handles.source,handles.pulse,...
        handles.Q_tsl,handles.f_tpp,handles.Q_tsl2,handles.dt,handles.maxt,handles.posx,...
        handles.posy,handles.posz,handles.gpumode,handles.N_abc,handles.b_abc,...
        handles.topoload,handles.topo,handles.X_topo,handles.Y_topo,...
        handles.savemode,handles.savetype,handles.skipt,handles.outpath);
elseif handles.stype==1 && handles.att_type==1 && handles.scheme==2
    % Point source, 4th order spatial derivatives, constant Q
    FD_sim_point_4(handles.Vp,handles.Vs,handles.rho,...
        handles.X,handles.Y,handles.Z,handles.source,handles.pulse,...
        handles.Q_tsl,handles.f_tpp,handles.dt,handles.maxt,handles.posx,...
        handles.posy,handles.posz,handles.gpumode,handles.N_abc,handles.b_abc,...
        handles.topoload,handles.topo,handles.X_topo,handles.Y_topo,...
        handles.savemode,handles.savetype,handles.skipt,handles.outpath);
elseif handles.stype==2 && handles.att_type==1 && handles.scheme==1
    % Fault source, 2nd order spatial derivatives, constant Q
    FD_sim_fault_2(handles.Vp,handles.Vs,handles.rho,...
        handles.X,handles.Y,handles.Z,handles.source,handles.pulse,...
        handles.Q_tsl,handles.f_tpp,handles.dt,handles.maxt,handles.posx,...
        handles.posy,handles.posz,handles.fault,handles.xrup,handles.zrup,handles.gpumode,...
        handles.N_abc,handles.b_abc,handles.topoload,handles.topo,handles.X_topo,handles.Y_topo,...
        handles.savemode,handles.savetype,handles.skipt,handles.outpath);
elseif handles.stype==2 && handles.att_type==2 && handles.scheme==1
    % Fault source, 2nd order spatial derivatives, Q based on relaxation functions
    FD_sim_fault_2_Q2(handles.Vp,handles.Vs,handles.rho,...
        handles.X,handles.Y,handles.Z,handles.source,handles.pulse,...
        handles.Q_tsl,handles.f_tpp,handles.Q_tsl2,handles.dt,handles.maxt,handles.posx,...
        handles.posy,handles.posz,handles.fault,handles.xrup,handles.zrup,handles.gpumode,...
        handles.N_abc,handles.b_abc,handles.topoload,handles.topo,handles.X_topo,handles.Y_topo,...
        handles.savemode,handles.savetype,handles.skipt,handles.outpath);
elseif handles.stype==2 && handles.att_type==2 && handles.scheme==2
    % Fault source, 4th order spatial derivatives, Q based on relaxation functions
    FD_sim_fault_4_Q2(handles.Vp,handles.Vs,handles.rho,...
        handles.X,handles.Y,handles.Z,handles.source,handles.pulse,...
        handles.Q_tsl,handles.f_tpp,handles.Q_tsl2,handles.dt,handles.maxt,handles.posx,...
        handles.posy,handles.posz,handles.fault,handles.xrup,handles.zrup,handles.gpumode,...
        handles.N_abc,handles.b_abc,handles.topoload,handles.topo,handles.X_topo,handles.Y_topo,...
        handles.savemode,handles.savetype,handles.skipt,handles.outpath);
elseif handles.stype==2 && handles.att_type==1 && handles.scheme==2
    % Fault source, 4th order spatial derivatives, constant Q
    FD_sim_fault_4(handles.Vp,handles.Vs,handles.rho,...
        handles.X,handles.Y,handles.Z,handles.source,handles.pulse,...
        handles.Q_tsl,handles.f_tpp,handles.dt,handles.maxt,handles.posx,...
        handles.posy,handles.posz,handles.fault,handles.xrup,handles.zrup,handles.gpumode,...
        handles.N_abc,handles.b_abc,handles.topoload,handles.topo,handles.X_topo,handles.Y_topo,...
        handles.savemode,handles.savetype,handles.skipt,handles.outpath);
end
