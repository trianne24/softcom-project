function varargout = eFSLab(varargin)
% EFSLAB M-file for eFSLab.fig
%      EFSLAB, by itself, creates a new EFSLAB or raises the existing
%      singleton*.
%
%      H = EFSLAB returns the handle to a new EFSLAB or the handle to
%      the existing singleton*.
%
%      EFSLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EFSLAB.M with the given input arguments.
%
%      EFSLAB('Property','Value',...) creates a new EFSLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Evolving_TS_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eFSLab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% The eFSLab has been developed at Center for Informatics and Systems of
% the University of Coimbra (CISUC). The first version was developed by
% José Victor Ramos and rhe Guide interface and second version of code was
% developed by Lara Aires, under the supervision of António Dourado. This work was partially suported by FCT project
% POSC/EIA/58162/2004 (CLASSE) envolving FEDER European Union funds.
%
% Copyright (C) 2008 António Dourado, Lara Aires, José Victor Ramos – FCTUC Coimbra, Portugal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To cite this work:
%
%           Report
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% Last Modified by GUIDE v2.5 06-Nov-2008 11:24:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eFSLab_OpeningFcn, ...
                   'gui_OutputFcn',  @eFSLab_OutputFcn, ...
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


% --- Executes just before eFSLab is made visible.


function eFSLab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eFSLab (see VARARGIN)

% Choose default command line output for eFSLab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes eFSLab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = eFSLab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in Browse_file_pushbutton.
function Browse_file_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_file_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 %[filename, path] = uigetfile({'*.*','All Files (*.*)';},'Select a file to load');
 
 %--- Read Excel and Text files 
 [filename, path] = uigetfile({'*.xls', 'Excel files(*.xls)';'*.txt','Text files (*.txt)';},'Select a file to load');

 if(filename~=0)
     set(handles.Select_file_listbox,'String',filename);
     handles.sourcefile.filename=filename;
     handles.sourcefile.path=path;
     set([handles.Select_file_listbox,
         handles.Input_columns_text,
         handles.Input_columns_edit,
         handles.Import_data_pushbutton],'enable','on');
 end
 
typ=xlsfinfo(handles.sourcefile.filename);
comp_str=strcmp(typ,'');

if comp_str==0
    set([handles.start_cell_text,
         handles.start_cell_edit,
         handles.end_cell_text,
         handles.end_cell_edit],'enable','on');
end
  
 
  % Update handles structure
  guidata(hObject, handles);    


% --- Executes on selection change in Select_file_listbox.
function Select_file_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Select_file_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Select_file_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Select_file_listbox


% --- Executes during object creation, after setting all properties.
function Select_file_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Select_file_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Input_columns_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Input_columns_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Input_columns_edit as text
%        str2double(get(hObject,'String')) returns contents of Input_columns_edit as a double


% --- Executes during object creation, after setting all properties.
function Input_columns_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_columns_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in Import_data_pushbutton.
function Import_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Import_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject,'Enable','Off')
pause(0.001);

%gets values in edit boxes and in list
    cell_s= get(handles.start_cell_edit,'String');
    cell_e= get(handles.end_cell_edit,'String');
    list_c=get(handles.Select_file_listbox,'String');
    val_c=get(handles.Select_file_listbox,'Value');
    input= get(handles.Input_columns_edit,'String');
    
% --- Read excel or text file    
typ=xlsfinfo(handles.sourcefile.filename);
comp_str=strcmp(typ,'');

if comp_str==1
    Init_data = load ([handles.sourcefile.path handles.sourcefile.filename]);
elseif comp_str==0
    Init_data=xlsread([handles.sourcefile.path handles.sourcefile.filename],strcat(cell_s,':',cell_e));
end

%////////////////// CASO DE FICHEIRO GENÉRICO  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
input=str2double(input);
N=input+1;
handles.N=N;
 
% Initialization of variables
[linhas colunas]=size(Init_data);
handles.linhas=linhas;
handles.colunas=colunas;
Dados = zeros(linhas,N);
Dados=[Init_data(:,1:input) Init_data(:,colunas)];
 
%Normalizacao dos dados de treino no intervalo [0 1]

for i=1:N
    min_tr(i) = min(Dados(:,i));
    max_tr(i) = max(Dados(:,i));
    
    if (max_tr(i)-min_tr(i))==0
        Dados(:,i)=Dados(:,i);
    else
        Dados(:,i) = (Dados(:,i)-min_tr(i))/(max_tr(i) - min_tr(i));
    end
end

handles.Data=Dados;
handles.input=input;

handles.Xin=Dados(:,1:(N-1));
handles.Xout=Dados(:,N);
guidata(hObject, handles);
set(hObject,'Enable','On')

%//////////////////////////////////// FIM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% %////////////////// CASO DOS DADOS Box-Jenkins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Dados = zeros(290, 3);
% 
% Dados = Init_data(:,[2,9,1]);
% N = 3;
% 
% % Normalizacao dos dados de treino no intervalo [0 1]
% for i=1:N
%      min_tr(i) = min(Dados(:,i));
%      max_tr(i) = max(Dados(:,i));
%      Dados(:,i) = (Dados(:,i)-min_tr(i))/(max_tr(i) - min_tr(i));
%  end
% 
% % Agregation of training and validation data sets
% % The data format will be [u(t-4), y(t-1); y(t)]
% Data = [Dados(:,2) Dados(:,1) Dados(:,3)];
% handles.Data=Data;
% handles.Xin=Data(:,1:(N-1));
% handles.Xout=Data(:,N);
% guidata(hObject, handles);
% set(hObject,'Enable','On')
% %//////////////////////////////////// FIM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Default of separator
[l c]=size(handles.Data);
sept=int2str(round(l*70/100));
set(handles.Separator_edit,'String',sept);


% --- Executes on selection change in Model_popupmenu.
function Model_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Model_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Model_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Model_popupmenu


% --- Executes during object creation, after setting all properties.
function Model_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Model_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in Parameter_Estimation_popupmenu.
function Parameter_Estimation_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Parameter_Estimation_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Parameter_Estimation_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Parameter_Estimation_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function Radii_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Radii_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Radii_edit as text
%        str2double(get(hObject,'String')) returns contents of Radii_edit as a double



% --- Executes during object creation, after setting all properties.
function Radii_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Radii_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Omega_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Omega_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Omega_edit as text
%        str2double(get(hObject,'String')) returns contents of Omega_edit as a double



% --- Executes during object creation, after setting all properties.
function Omega_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Omega_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Delta_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Delta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Delta_edit as text
%        str2double(get(hObject,'String')) returns contents of Delta_edit as a double


% --- Executes during object creation, after setting all properties.
function Delta_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Delta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Epson_min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Epson_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Epson_min_edit as text
%        str2double(get(hObject,'String')) returns contents of Epson_min_edit as a double


% --- Executes during object creation, after setting all properties.
function Epson_min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Epson_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Epson_max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Epson_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Epson_max_edit as text
%        str2double(get(hObject,'String')) returns contents of Epson_max_edit as a double


% --- Executes during object creation, after setting all properties.
function Epson_max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Epson_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Separator_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Separator_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Separator_edit as text
%        str2double(get(hObject,'String')) returns contents of Separator_edit as a double


% --- Executes during object creation, after setting all properties.
function Separator_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Separator_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end






% --- Executes on selection change in Conditions_Substitution_Rules_listbox.
function Conditions_Substitution_Rules_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Conditions_Substitution_Rules_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Conditions_Substitution_Rules_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Conditions_Substitution_Rules_listbox


% --- Executes during object creation, after setting all properties.
function Conditions_Substitution_Rules_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conditions_Substitution_Rules_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in Conditions_Create_Rules_listbox.
function Conditions_Create_Rules_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Conditions_Create_Rules_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Conditions_Create_Rules_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Conditions_Create_Rules_listbox



% --- Executes during object creation, after setting all properties.
function Conditions_Create_Rules_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conditions_Create_Rules_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Reference_2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Reference_2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Reference_2_edit as text
%        str2double(get(hObject,'String')) returns contents of Reference_2_edit as a double



% --- Executes during object creation, after setting all properties.
function Reference_2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Reference_2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Reference_1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Reference_1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Reference_1_edit as text
%        str2double(get(hObject,'String')) returns contents of Reference_1_edit as a double


% --- Executes during object creation, after setting all properties.
function Reference_1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Reference_1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function My_creation_rule_condition_edit_Callback(hObject, eventdata, handles)
% hObject    handle to My_creation_rule_condition_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of My_creation_rule_condition_edit as text
%        str2double(get(hObject,'String')) returns contents of My_creation_rule_condition_edit as a double



% --- Executes during object creation, after setting all properties.
function My_creation_rule_condition_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to My_creation_rule_condition_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function My_Substitution_Rules_Condition_edit_Callback(hObject, eventdata, handles)
% hObject    handle to My_Substitution_Rules_Condition_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of My_Substitution_Rules_Condition_edit as text
%        str2double(get(hObject,'String')) returns contents of My_Substitution_Rules_Condition_edit as a double



% --- Executes during object creation, after setting all properties.
function My_Substitution_Rules_Condition_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to My_Substitution_Rules_Condition_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ok_pushbutton.
function ok_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ok_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Conditions_Substitution_Rules_listbox_str=get(handles.Conditions_Substitution_Rules_listbox,'String');
Conditions_Substitution_Rules_listbox_val=get(handles.Conditions_Substitution_Rules_listbox,'Value');

switch Conditions_Substitution_Rules_listbox_str{Conditions_Substitution_Rules_listbox_val};
    case '1. ((Potential(k) > Potential_Ref(k))  & (Distance_Min / radii < Reference_1))  |  ((Potential(k) < Potential_Min(k))  &  (Distance_Min / radii < Reference_2))  '
        available_create_conditions=cellstr(strvcat('1. (Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k))','2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);
    case '2. (Potential(k) > Epsilon_down(k)) & (Distance_Min / radii + Potential_Ref(k) / Potential(k)) < Reference_1'
        available_create_conditions=cellstr(strvcat('2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);
    case '3. (Potential(k) > Potential_Mean(k)) & (Distance_Min / radii + Potential_Ref(k) / Potential(k)) < Reference_1'
        available_create_conditions=cellstr(strvcat('2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);
    case '4. (Potential(k) > Potential_Ref(k)) & (Distance_Min / radii + Potential_Ref(k) / Potential(k)) < Reference_1'
        available_create_conditions=cellstr(strvcat('2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);
    case '5. (Potential(k) > Potential_Ref(k)) & (Distance_Min / radii) < Reference_1'
        available_create_conditions=cellstr(strvcat('1. (Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k))','2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);        
    case '6. ((Potential(k) > Potential_Ref(k))  &  (Distance_Min / radii < Reference_1))  | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)) & (Distance_Min / radii < Reference_1))'
        available_create_conditions=cellstr(strvcat('1. (Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k))','2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions); 
    case '7. ((Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k)))  &  (Distance_Min / radii < Reference_1)' 
        available_create_conditions=cellstr(strvcat('1. (Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k))','2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);    
    case '8. (Potential(k) > Potential_Ref(k)) &  (Distance_Min / radii  +  Potential_Ref(k)  /  Potential(k)) < Reference_1'
        available_create_conditions=cellstr(strvcat('2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);        
    case '9. (Potential(k) > Potential_Ref(k))  &  (Distance_Min / radii  +  Potential_Mean(k) /  Potential(k)) < Reference_1'
        available_create_conditions=cellstr(strvcat('2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);       
    case '10. ((Potential(k) > Potential_Mean(k))  |  ((Potential(k) > Epsilon_down(k))  &  (Potential(k) < Epsilon_up(k))))  &  (Distance_Min / radii < Reference_1)'
        available_create_conditions=cellstr(strvcat('2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))','4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);  
    case '11. ((Potential(k) > Potential_Ref(k))  &  (Distance_Min / radii < Reference_1)) | ((Potential(k) < Potential_Min(k)) & (Distance_Min / radii + Potential(k) / Potential_Min(k)) >= Reference_2) '
        available_create_conditions=cellstr(strvcat('1. (Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k))','2. (Potential(k) > Potential_Mean(k))','3. (Potential(k) > Potential_Ref(k))'));
        set(handles.Conditions_Create_Rules_listbox,'String',available_create_conditions);        
end


% --- Executes on button press in Create_Model_pushbutton.
function Create_Model_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Create_Model_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


model_str=get(handles.Model_popupmenu,'String');
model_val=get(handles.Model_popupmenu,'Value');

switch model_str{model_val};
    case 'Takagi-Sugeno 0 order'
        handles.current_model=0;
    case 'Takagi-Sugeno 1st order'
        handles.current_model=1;
end

parameter_est_str=get(handles.Parameter_Estimation_popupmenu,'String');
parameter_est_val=get(handles.Parameter_Estimation_popupmenu,'Value');

switch parameter_est_str{parameter_est_val};
    case 'RLS - Global Estimation'
        handles.current_estimation=1;
    case 'wRLS - Local Estimation'
        handles.current_estimation=2;
end

handles.radii=str2double(get(handles.Radii_edit,'String'));
handles.Omega=str2double(get(handles.Omega_edit,'String'));
handles.delta=str2double(get(handles.Delta_edit,'String'));
handles.epson_down=str2double(get(handles.Epson_min_edit,'String'));
handles.epson_up=str2double(get(handles.Epson_max_edit,'String'));
handles.separador=str2double(get(handles.Separator_edit,'String'));


Conditions_Substitution_Rules_listbox_str=get(handles.Conditions_Substitution_Rules_listbox,'String');
Conditions_Substitution_Rules_listbox_val=get(handles.Conditions_Substitution_Rules_listbox,'Value');

switch Conditions_Substitution_Rules_listbox_str{Conditions_Substitution_Rules_listbox_val};
    case '1. ((Potential(k) > Potential_Ref(k))  & (Distance_Min / radii < Reference_1))  |  ((Potential(k) < Potential_Min(k))  &  (Distance_Min / radii < Reference_2))  '
        handles.current_substitution_condition=' ((Potential(k) > Potential_Ref(k))  & (Distance_Min / radii < Reference_1))  |  ((Potential(k) < Potential_Min(k))  &  (Distance_Min / radii < Reference_2))  ';
    case '2. (Potential(k) > Epsilon_down(k)) & (Distance_Min / radii + Potential_Ref(k) / Potential(k)) < Reference_1'
        handles.current_substitution_condition=' (Potential(k) > Epsilon_down(k)) & (Distance_Min / radii + Potential_Ref(k) / Potential(k)) < Reference_1';
    case '3. (Potential(k) > Potential_Mean(k)) & (Distance_Min / radii + Potential_Ref(k) / Potential(k)) < Reference_1'
        handles.current_substitution_condition=' (Potential(k) > Potential_Mean(k)) & (Distance_Min / radii + Potential_Ref(k) / Potential(k)) < Reference_1';
    case '4. (Potential(k) > Potential_Ref(k)) & (Distance_Min / radii + Potential_Ref(k) / Potential(k)) < Reference_1'
        handles.current_substitution_condition=' (Potential(k) > Potential_Ref(k)) & (Distance_Min / radii + Potential_Ref(k) / Potential(k)) < Reference_1';
    case '5. (Potential(k) > Potential_Ref(k)) & (Distance_Min / radii) < Reference_1'
        handles.current_substitution_condition=' (Potential(k) > Potential_Ref(k)) & (Distance_Min / radii) < Reference_1';
    case '6. ((Potential(k) > Potential_Ref(k))  &  (Distance_Min / radii < Reference_1))  | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)) & (Distance_Min / radii < Reference_1))'
        handles.current_substitution_condition=' ((Potential(k) > Potential_Ref(k))  &  (Distance_Min / radii < Reference_1))  | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)) & (Distance_Min / radii < Reference_1))';
    case '7. ((Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k)))  &  (Distance_Min / radii < Reference_1)' 
        handles.current_substitution_condition=' ((Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k)))  &  (Distance_Min / radii < Reference_1)' ;
    case '8. (Potential(k) > Potential_Ref(k)) &  (Distance_Min / radii  +  Potential_Ref(k)  /  Potential(k)) < Reference_1'
        handles.current_substitution_condition=' (Potential(k) > Potential_Ref(k)) &  (Distance_Min / radii  +  Potential_Ref(k)  /  Potential(k)) < Reference_1';
    case '9. (Potential(k) > Potential_Ref(k))  &  (Distance_Min / radii  +  Potential_Mean(k) /  Potential(k)) < Reference_1'
        handles.current_substitution_condition=' (Potential(k) > Potential_Ref(k))  &  (Distance_Min / radii  +  Potential_Mean(k) /  Potential(k)) < Reference_1';
    case '10. ((Potential(k) > Potential_Mean(k))  |  ((Potential(k) > Epsilon_down(k))  &  (Potential(k) < Epsilon_up(k))))  &  (Distance_Min / radii < Reference_1)'
        handles.current_substitution_condition=' ((Potential(k) > Potential_Mean(k))  |  ((Potential(k) > Epsilon_down(k))  &  (Potential(k) < Epsilon_up(k))))  &  (Distance_Min / radii < Reference_1)';
    case '11. ((Potential(k) > Potential_Ref(k))  &  (Distance_Min / radii < Reference_1)) | ((Potential(k) < Potential_Min(k)) & (Distance_Min / radii + Potential(k) / Potential_Min(k)) >= Reference_2) '
        handles.current_substitution_condition=' ((Potential(k) > Potential_Ref(k))  &  (Distance_Min / radii < Reference_1)) | ((Potential(k) < Potential_Min(k)) & (Distance_Min / radii + Potential(k) / Potential_Min(k)) >= Reference_2) ';
end

Conditions_Create_Rules_listbox_str2=get(handles.Conditions_Create_Rules_listbox,'String');
Conditions_Create_Rules_listbox_val2=get(handles.Conditions_Create_Rules_listbox,'Value');

switch Conditions_Create_Rules_listbox_str2{Conditions_Create_Rules_listbox_val2};
    case '1. (Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k))'
        handles.current_create_rules_condition=' (Potential(k) > Potential_Ref(k)) | (Potential(k) < Potential_Min(k))';
    case '2. (Potential(k) > Potential_Mean(k))'
        handles.current_create_rules_condition=' (Potential(k) > Potential_Mean(k))';
    case '3. (Potential(k) > Potential_Ref(k))'
        handles.current_create_rules_condition=' (Potential(k) > Potential_Ref(k))';
    case '4. (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))'
        handles.current_create_rules_condition=' (Potential(k) > Potential_Ref(k)) | ((Potential(k) > Epsilon_down(k)) & (Potential(k) < Epsilon_up(k)))';
end

handles.My_creation_rule_condition=get(handles.My_creation_rule_condition_edit,'String');
handles.My_Substitution_Rules_Condition=get(handles.My_Substitution_Rules_Condition_edit,'String');

handles.Reference_1=str2double(get(handles.Reference_1_edit,'String'));
handles.Reference_2=str2double(get(handles.Reference_2_edit,'String'));


if length(char(handles.My_Substitution_Rules_Condition))>0
    handles.current_substitution_condition=handles.My_Substitution_Rules_Condition;
end

if length(char(handles.My_creation_rule_condition))>0
    handles.current_create_rules_condition=handles.My_creation_rule_condition;
end

if (get(handles.View_Graphics_checkbox,'Value') == get(handles.View_Graphics_checkbox,'Max'))
	% Checkbox is checked-take approriate action
    graphics=1;
    handles.graphics=graphics;
else
	graphics=0;
    handles.graphics=graphics;
end

[centers, sigmas, conseq, R] = evolving_TS(handles.input, handles.Data, handles.radii, handles.Omega, handles.delta, handles.epson_down, handles.epson_up, handles.separador,handles.Reference_1, handles.Reference_2, handles.current_substitution_condition, handles.current_create_rules_condition, handles.current_model, handles.current_estimation,  handles.graphics);

 for i=1:R
     sigmas_f(i,:)=sigmas;
 end

handles.centers=centers;
handles.sigmas=sigmas;
handles.conseq=conseq;
handles.conseq_t=(conseq)';
handles.sigmas_f=sigmas_f;
handles.R=R;

fismat= modelo_rsc(handles.Xin, handles.Xout, handles.centers, handles.sigmas, handles.conseq, handles.current_model);
handles.Sug_fismat=fismat;

set(handles.Creat_Results_Table_pushbutton,'Enable','on');
set(handles.Transform_Sugeno_into_Mamdani_pushbutton,'Enable','on');
set(handles.Open_in_Fuzzy_Toolbox_pushbutton,'Enable','on');
set(handles.Reset_pushbutton,'Enable','on');

guidata(hObject,handles);


function Apresentacao_ficheiro_excel_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Apresentacao_ficheiro_excel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Apresentacao_ficheiro_excel_edit as text
%        str2double(get(hObject,'String')) returns contents of Apresentacao_ficheiro_excel_edit as a double


% --- Executes during object creation, after setting all properties.
function Apresentacao_ficheiro_excel_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Apresentacao_ficheiro_excel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Creat_Results_Table_pushbutton.
function Creat_Results_Table_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Creat_Results_Table_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handle_uitable = uitable ('Parent',handles.Rules_Antecedents_and_Consequents_uipanel,'Tag','Results_uitable','Units','normalized');

handles.handle_uitable=handle_uitable;
rules=[];

for i=1:(handles.R)
    num=int2str(i);
    str=strcat('Rule',num);
    rules=strvcat(rules,str);
end

row_name=rules;

if handles.current_model == 0
    column_name= strvcat('AntC1','AntC2','AntSig1','AntSig2','ConseqC0');
    dados_f=horzcat(handles.centers, handles.sigmas_f, handles.conseq_t);
elseif handles.current_model == 1
    column_name= strvcat('AntC1','AntC2','AntSig1','AntSig2','ConseqC0','ConseqC1(u1)','ConseqC2(u2)');
    dados_f=horzcat(handles.centers, handles.sigmas_f, handles.conseq);
end

set(handle_uitable,'RowName',row_name,'ColumnName',column_name,'Data',dados_f);

pos=get(handles.Rules_Antecedents_and_Consequents_uipanel,'Position');
pos(1)=0.025;
pos(2)=0.1;
pos(3)=pos(3)-0.005;
pos(4)=pos(4)+0.05;
set(handle_uitable,'Position',(2*pos));

guidata(hObject,handles);


function start_cell_edit_Callback(hObject, eventdata, handles)
% hObject    handle to start_cell_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_cell_edit as text
%        str2double(get(hObject,'String')) returns contents of start_cell_edit as a double


% --- Executes during object creation, after setting all properties.
function start_cell_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_cell_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function end_cell_edit_Callback(hObject, eventdata, handles)
% hObject    handle to end_cell_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_cell_edit as text
%        str2double(get(hObject,'String')) returns contents of end_cell_edit as a double


% --- Executes during object creation, after setting all properties.
function end_cell_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_cell_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Open_in_Fuzzy_Toolbox_pushbutton.
function Open_in_Fuzzy_Toolbox_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Open_in_Fuzzy_Toolbox_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fuzzy(handles.Sug_fismat);
guidata(hObject, handles);

% --- Executes on button press in Transform_Sugeno_into_Mamdani_pushbutton.
function Transform_Sugeno_into_Mamdani_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Transform_Sugeno_into_Mamdani_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.current_model==1
    warndlg('Mamdani precision would be worst than TSK. You will transform first order TSK into Mamdani, centering the consequent in the corresponding TSK independent consequent.',' Warning ');
    uiwait
end

SugenoToMamdani('eFSLab', handles.figure1);
%guidata(hObject, handles);



% --- Executes on button press in Reset_pushbutton.
function Reset_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Data
% set(handles.Select_file_listbox,'String','');
% set(handles.start_cell_text,'Enable','off');
% set(handles.start_cell_edit,'Enable','off');
% set(handles.start_cell_edit,'String','C2');
% set(handles.end_cell_text,'Enable','off');
% set(handles.end_cell_edit,'Enable','off');
% set(handles.end_cell_edit,'String','AW2884');
% set(handles.Input_columns_edit,'String','2');
% 
% % Model & Estimation
% set(handles.Model_popupmenu,'Value',1.0);
% set(handles.Parameter_Estimation_popupmenu,'Value',1.0);
% 
% % Other Parameters Values
% set(handles.Radii_edit,'String','0.3');
% set(handles.Omega_edit,'String','750');
% set(handles.Delta_edit,'String','0.1');
% set(handles.Epson_min_edit,'String','0.5');
% set(handles.Epson_max_edit,'String','0.75');
% set(handles.Separator_edit,'String','');
% 
% % Substitution and Creation of Rules
% set(handles.Conditions_Substitution_Rules_listbox,'Value',1.0);
% set(handles.Reference_1_edit,'String','0.5');
% set(handles.Reference_2_edit,'String','0.85');
% 
% set(handles.Conditions_Create_Rules_listbox,'Value',1.0);
% 
% % My Conditions
% set(handles.My_Substitution_Rules_Condition_edit,'String','');
% set(handles.My_creation_rule_condition_edit,'String','');
% 
% % Check Box
% set(handles.View_Graphics_checkbox,'Value',0)
% 
% % Buttons
% set(handles.Creat_Results_Table_pushbutton,'Enable','off');
% set(handles.Open_in_Fuzzy_Toolbox_pushbutton,'Enable','off');
% set(handles.Transform_Sugeno_into_Mamdani_pushbutton,'Enable','off');
% 
% % Table
% Hmatch_t = findobj('Tag','Results_uitable');
% emp_Hmatch_t=isempty(Hmatch_t);
% if emp_Hmatch_t ==0
%     delete(handles.handle_uitable);
% end
% 
% % Reset
% set(handles.Reset_pushbutton,'Enable','off');
% 
% % Handles
% Hmatch_f = findobj('Name','SugenoToMamdani');
% emp_Hmatch_f=isempty(Hmatch_f);
% if emp_Hmatch_f == 0
%     close SugenoToMamdani;
% end

close all;
clear all;
eFSLab;



% --- Executes on button press in View_Graphics_checkbox.
function View_Graphics_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to View_Graphics_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of View_Graphics_checkbox

% if (get(hObject,'Value') == get(hObject,'Max'))
% 	% Checkbox is checked-take approriate action
%     graphics=1;
%     handles.graphics=graphics;
% else
% 	graphics=0;
%     handles.graphics=graphics;
% end
% 
% guidata(hObject, handles);



