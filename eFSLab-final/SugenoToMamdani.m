function varargout = SugenoToMamdani(varargin)
% SUGENOTOMAMDANI M-file for SugenoToMamdani.fig
%      SUGENOTOMAMDANI, by itself, creates a new SUGENOTOMAMDANI or raises the existing
%      singleton*.
%
%      H = SUGENOTOMAMDANI returns the handle to a new SUGENOTOMAMDANI or the handle to
%      the existing singleton*.
%
%      SUGENOTOMAMDANI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUGENOTOMAMDANI.M with the given input arguments.
%
%      SUGENOTOMAMDANI('Property','Value',...) creates a new SUGENOTOMAMDANI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SugenoToMamdani_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SugenoToMamdani_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SugenoToMamdani

% Last Modified by GUIDE v2.5 29-Oct-2008 17:03:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SugenoToMamdani_OpeningFcn, ...
                   'gui_OutputFcn',  @SugenoToMamdani_OutputFcn, ...
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


% --- Executes just before SugenoToMamdani is made visible.
function SugenoToMamdani_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SugenoToMamdani (see VARARGIN)

dontOpen = false;
mainGuiInput = find(strcmp(varargin, 'eFSLab'));
if (isempty(mainGuiInput)) || (length(varargin) <= mainGuiInput) || (~ishandle(varargin{mainGuiInput+1}))
    dontOpen = true;
else
    
    % Remember the handle
    handles.eFSLabMain = varargin{mainGuiInput+1};
end

% Update handles structure
guidata(hObject, handles);
uiwait(hObject);

% UIWAIT makes SugenoToMamdani wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%--- Outputs from this function are returned to the command line.
function varargout = SugenoToMamdani_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;

varargout{1} = [];
delete(hObject);



function Partition_Number_gbellmf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Partition_Number_gbellmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Partition_Number_gbellmf_edit as text
%        str2double(get(hObject,'String')) returns contents of Partition_Number_gbellmf_edit as a double


% --- Executes during object creation, after setting all properties.
function Partition_Number_gbellmf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Partition_Number_gbellmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Slope_gbellmf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Slope_gbellmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Slope_gbellmf_edit as text
%        str2double(get(hObject,'String')) returns contents of Slope_gbellmf_edit as a double


% --- Executes during object creation, after setting all properties.
function Slope_gbellmf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slope_gbellmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Partition_Number_trimf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Partition_Number_trimf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Partition_Number_trimf_edit as text
%        str2double(get(hObject,'String')) returns contents of Partition_Number_trimf_edit as a double


% --- Executes during object creation, after setting all properties.
function Partition_Number_trimf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Partition_Number_trimf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Left_Slope_dsigmf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Left_Slope_dsigmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Left_Slope_dsigmf_edit as text
%        str2double(get(hObject,'String')) returns contents of Left_Slope_dsigmf_edit as a double


% --- Executes during object creation, after setting all properties.
function Left_Slope_dsigmf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Left_Slope_dsigmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Partition_Number_dsigmf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Partition_Number_dsigmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Partition_Number_dsigmf_edit as text
%        str2double(get(hObject,'String')) returns contents of Partition_Number_dsigmf_edit as a double


% --- Executes during object creation, after setting all properties.
function Partition_Number_dsigmf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Partition_Number_dsigmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Right_Slope_dsigmf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Right_Slope_dsigmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Right_Slope_dsigmf_edit as text
%        str2double(get(hObject,'String')) returns contents of Right_Slope_dsigmf_edit as a double


% --- Executes during object creation, after setting all properties.
function Right_Slope_dsigmf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Right_Slope_dsigmf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Preview_pushbutton.
function Preview_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Preview_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
option_sel=get(get(handles.Membership_Function_Type_uipanel,'SelectedObject'),'String');

main=handles.eFSLabMain;
mainHandles = guidata(main);

switch option_sel
    case 'trimf'
        part_numb=str2double(get(handles.Partition_Number_trimf_edit,'String'));
        conseqs_main=mainHandles.conseq;
                
        x=0:0.01:1;
        aperture=(1/(part_numb-1));
        params=[(conseqs_main(3)-aperture) conseqs_main(3) (conseqs_main(3)+aperture)];
        
        y=trimf(x,params);
        plot(handles.Preview_axes,x,y);
        
    case 'gaussmf'
        conseqs_main=mainHandles.conseq;
        sigmas_main=mainHandles.sigmas;
        
        params=[sigmas_main(1) conseqs_main(3)];
        x=0:0.01:1;
        y=gaussmf(x,params);
        
        plot(handles.Preview_axes,x,y);

    case 'gbellmf'
        part_numb=str2double(get(handles.Partition_Number_gbellmf_edit,'String'));
        slope=str2double(get(handles.Slope_gbellmf_edit,'String'));
        conseqs_main=mainHandles.conseq;
        aperture=1/part_numb;
        
        params=[aperture slope conseqs_main(3)];
        x=0:0.01:1;
        y=gbellmf(x,params);
        
        plot(handles.Preview_axes,x,y);

    case 'dsigmf'
        left_slope=str2double(get(handles.Left_Slope_dsigmf_edit,'String'));
        part_numb=str2double(get(handles.Partition_Number_dsigmf_edit,'String'));
        right_slope=str2double(get(handles.Right_Slope_dsigmf_edit,'String')); 
        conseqs_main=mainHandles.conseq;
        
        left_aperture=conseqs_main(3)-(1*2/part_numb);
        right_aperture=conseqs_main(3)+(1*2/part_numb);
        
        params=[left_slope left_aperture right_slope right_aperture];
        x=0:0.01:1;
        y=dsigmf(x,params);
        
        plot(handles.Preview_axes,x,y);
           
end 
guidata(hObject, handles);


% --- Executes on button press in Transform_pushbutton.
function Transform_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Transform_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mf_type=get(get(handles.Membership_Function_Type_uipanel,'SelectedObject'),'String');

main=handles.eFSLabMain;
mainHandles = guidata(main);
conseqs_main=mainHandles.conseq;
sigmas_main=mainHandles.sigmas;
fismat_main=mainHandles.Sug_fismat;

part_numb_trimf=str2double(get(handles.Partition_Number_trimf_edit,'String'));

part_numb_gbellmf=str2double(get(handles.Partition_Number_gbellmf_edit,'String'));
slope_gbellmf=str2double(get(handles.Slope_gbellmf_edit,'String'));

left_slope_dsigmf=str2double(get(handles.Left_Slope_dsigmf_edit,'String'));
part_numb_dsigmf=str2double(get(handles.Partition_Number_dsigmf_edit,'String'));
right_slope_dsigmf=str2double(get(handles.Right_Slope_dsigmf_edit,'String')); 

out_fismat = sug2mam(fismat_main,sigmas_main,conseqs_main,mf_type,part_numb_trimf, part_numb_gbellmf, slope_gbellmf, left_slope_dsigmf, part_numb_dsigmf, right_slope_dsigmf);
handles.out_fismat=out_fismat;

set(handles.Open_Fuzzy_Toolbox_pushbutton,'Enable','on');

guidata(hObject, handles);


% --- Executes on button press in Open_Fuzzy_Toolbox_pushbutton.
function Open_Fuzzy_Toolbox_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Open_Fuzzy_Toolbox_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%fismat=handles.Mam_fismat;
fuzzy(handles.out_fismat);

guidata(hObject, handles);


% --- Executes when selected object is changed in Membership_Function_Type_uipanel.
function Membership_Function_Type_uipanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Membership_Function_Type_uipanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
option=get(get(handles.Membership_Function_Type_uipanel,'SelectedObject'),'String');
switch option
    case 'trimf'
        set(handles.Partition_Number_trimf_text,'Enable','on');
        set(handles.Partition_Number_trimf_edit,'Enable','on');
    case 'gaussmf'
        set(handles.No_parameters_defined_text,'Enable','on');
    case 'gbellmf'
        set(handles.Partition_Number_gbellmf_text,'Enable','on');
        set(handles.Partition_Number_gbellmf_edit,'Enable','on');
        
        set(handles.Slope_gbellmf_text,'Enable','on');
        set(handles.Slope_gbellmf_edit,'Enable','on');
    case 'dsigmf'
        set(handles.Left_Slope_dsigmf_text,'Enable','on');
        set(handles.Left_Slope_dsigmf_edit,'Enable','on');
        
        set(handles.Partition_Number_dsigmf_text,'Enable','on');
        set(handles.Partition_Number_dsigmf_edit,'Enable','on');
        
        set(handles.Right_Slope_dsigmf_text,'Enable','on');
        set(handles.Right_Slope_dsigmf_edit,'Enable','on');    
end 

oldValue=get(eventdata.OldValue,'String');
newValue=get(eventdata.NewValue,'String');

if strcmp(oldValue,newValue)==0
    if strcmp(oldValue,'trimf')
            set(handles.Partition_Number_trimf_text,'Enable','off');
            set(handles.Partition_Number_trimf_edit,'Enable','off');
            
    elseif strcmp(oldValue,'gaussmf')
            set(handles.No_parameters_defined_text,'Enable','off');
            
    elseif strcmp(oldValue,'gbellmf')            
            set(handles.Partition_Number_gbellmf_text,'Enable','off');
            set(handles.Partition_Number_gbellmf_edit,'Enable','off');
            
            set(handles.Slope_gbellmf_text,'Enable','off');
            set(handles.Slope_gbellmf_edit,'Enable','off');
            
    elseif strcmp(oldValue,'dsigmf')   
            set(handles.Left_Slope_dsigmf_text,'Enable','off');
            set(handles.Left_Slope_dsigmf_edit,'Enable','off');
            
            set(handles.Partition_Number_dsigmf_text,'Enable','off');
            set(handles.Partition_Number_dsigmf_edit,'Enable','off');
            
            set(handles.Right_Slope_dsigmf_text,'Enable','off');
            set(handles.Right_Slope_dsigmf_edit,'Enable','off'); 
    end
end
guidata(hObject, handles);

    
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(hObject);


