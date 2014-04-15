function varargout = removePLI_gui(varargin)
% Graphical User Interface
% for Power Line Interference Cancellation for multi channel data
%   This is an implementation of the proposed algorithm in,
%   M. R. Keshtkaran and Z. Yang, "A fast, robust algorithm for power line 
%   interference cancellation in neural recording," J. Neural Eng., vol. 11,
%   no. 2, p. 026017, Apr. 2014.
%  for more information type: help remove_PLI_multichan

%   Licence:
%   Downloaded from: https://github.com/mrezak/removePLI
%	Author: Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
%   Copyright (c) 2013, Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
%   All rights reserved.
%	This program is provided "AS IS" for non-commercial, educational 
%	and reseach purpose only. Any commercial use, of any kind, of 
%	this program is prohibited. The Copyright notice should remain intact.
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @removePLI_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @removePLI_gui_OutputFcn, ...
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

function runRemovePLI(handles)
    contents = cellstr(get(handles.popInputVar,'String')); % returns popInputVar contents as cell array
    InVarN = contents{get(handles.popInputVar,'Value')}; %returns selected item from popInputVar
    OutVar = get(handles.edtOutputVar,'String');
    X = evalin('base',InVarN);
    szX = size(squeeze(X));
    if length(szX)~=2 || szX(2)<2 
        msgbox('Inappropriate input. The input must be a M*N matrix where M is the number of channels and N is the number of samples.',...
        'Error: Inappropriate input', 'error', 'modal')
        return
    end
    pChanFreq = get(handles.popChanFreq,'Value')-1;
    pSR = str2double(get(handles.edtSRate,'String'));
    if ~(pSR > 1)
        msgbox('Please input the correct sampling rate.','Error: Samling rate','error', 'modal')
        return
    end
    switch get(handles.popACFreq,'Value')
        case 1
            f_ac = [];
        case 2
            f_ac = 50;
        case 3 
            f_ac = 60;
    end
    contents = cellstr(get(handles.popHarm,'String')); % returns popInputVar contents as cell array
    M = str2double(contents{get(handles.popHarm,'Value')}); %returns selected item from popInputVar
    B0 = str2double(get(handles.edtB0,'String'));
    Binf = str2double(get(handles.edtBinf,'String'));
    Bst = str2double(get(handles.edtBst,'String'));
    P0 = str2double(get(handles.edtP0,'String'));
    Pinf = str2double(get(handles.edtPinf,'String'));
    Pst = str2double(get(handles.edtPst,'String'));
    W = str2double(get(handles.edtW,'String'));
    
    % Performing Interference Cancellation
    Y = removePLI_multichan(X, pSR, M,[B0 Binf Bst], [P0 Pinf Pst], W, f_ac, pChanFreq);
    assignin('base',OutVar,Y);
    % Plotting the output
    set(handles.popPlotChan,'Enable','on')
    set(handles.edtSkip,'Enable','on')
    set(handles.btnPlot,'Enable', 'on');
    DoPlots(handles)
    
    
function DoPlots(handles)
    contents = cellstr(get(handles.popInputVar,'String')); % returns popInputVar contents as cell array
    InVarN = contents{get(handles.popInputVar,'Value')}; %returns selected item from popInputVar
    X = evalin('base',InVarN);
    OutVarN = get(handles.edtOutputVar,'String');
    Y = evalin('base',OutVarN);
    pSR = str2double(get(handles.edtSRate,'String'));
    contents = cellstr(get(handles.popPlotChan,'String')); % returns popInputVar contents as cell array
    ChanNum = str2double(contents{get(handles.popPlotChan,'Value')});
    SkipInit = str2double(get(handles.edtSkip,'String'));
    st = min(pSR*SkipInit,(size(X,2)-30*pSR));
    ed = st + 30*pSR;
    idx = st:ed; %select 30 seconds of the data for plotting
    idx = idx(idx>0);
    axes(handles.axesBefore);
    pwelch(X(ChanNum,idx),[],[],[],pSR);
    s = sprintf('PSD of Channel %d Before Inteference Cancellation (30 sec Segment)',ChanNum);
    title(s);
    axes(handles.axesAfter);
    pwelch(Y(ChanNum,idx),[],[],[],pSR);
    s = sprintf('PSD of Channel %d After Inteference Cancellation (30 sec Segment)',ChanNum);
    title(s);
    

function removePLI_gui_InitFields(handles)
    varnames = evalin('base','who');
    set(handles.popInputVar,'String',[varnames]);
    contents = cellstr(get(handles.popInputVar,'String')); % returns popInputVar contents as cell array
    InVarN = contents{get(handles.popInputVar,'Value')}; %returns selected item from popInputVar
    set(handles.edtOutputVar,'String',[InVarN '_clean'])    
    InVar = evalin('base',InVarN);
    szX = size(InVar);
    set(handles.popPlotChan,'Enable', 'off');
    set(handles.edtSkip,'Enable', 'off');
    set(handles.btnPlot,'Enable', 'off');
    if length(szX)~=2 || szX(2)<2 
       set(handles.popChanFreq,'Value',1)
       set(handles.popPlotChan,'Value',1)
       set(handles.popHarm,'Value',1)
       set(handles.txtChan,'String', '---');
       set(handles.txtSamp,'String', '---');
       set(handles.popChanFreq,'String', '---');        
       set(handles.popPlotChan,'String', '---');        
       set(handles.txtChan,'Enable', 'off');
       set(handles.txtSamp,'Enable', 'off');
       set(handles.popChanFreq,'Enable', 'off');        
       set(handles.popChanFreq,'Enable', 'off');
       set(handles.edtSRate,'Enable', 'off');
       set(handles.popHarm,'Enable', 'off');
       set(handles.popACFreq,'Enable', 'off');
       set(handles.edtOutputVar,'Enable', 'off');
         return
    end
       set(handles.txtChan,'Enable', 'on');
       set(handles.txtSamp,'Enable', 'on');
       set(handles.popChanFreq,'Enable', 'on');        
       set(handles.popChanFreq,'Enable', 'on');
       set(handles.edtSRate,'Enable', 'on');
       set(handles.popHarm,'Enable', 'on');
       set(handles.popACFreq,'Enable', 'on');
       set(handles.edtOutputVar,'Enable', 'on');
       set(handles.txtChan,'String', num2str(size(InVar,1)));
       set(handles.txtSamp,'String', num2str(size(InVar,2)));
       set(handles.popChanFreq,'String', [num2cell(1:size(InVar,1)),{'Individual'}]);
       set(handles.popPlotChan,'String', num2cell(1:size(InVar,1)));
        linkaxes([handles.axesBefore handles.axesAfter],'xy');
        set(handles.popChanFreq,'Value',1)
        set(handles.popPlotChan,'Value',1)
        set(handles.edtSRate,'String', '');
    if strcmp(InVarN, 'gui_example'),
        set(handles.edtSRate,'String', '500');
        edtSRate_Callback(handles.edtSRate, [], handles)
        set(handles.popHarm,'Value', 3);
    end
        
        InitPSD(handles)

        
function InitPSD(handles)
    contents = cellstr(get(handles.popInputVar,'String')); % returns popInputVar contents as cell array
    InVarN = contents{get(handles.popInputVar,'Value')}; %returns selected item from popInputVar
    contents = cellstr(get(handles.popChanFreq,'String')); % returns popInputVar contents as cell array
    ChanFreq = str2double(contents{get(handles.popChanFreq,'Value')}); %returns selected item from popInputVar
    X = evalin('base',InVarN);
    axes(handles.axesBefore);
    pSR = str2double(get(handles.edtSRate,'String'));
    cla(handles.axesAfter,'reset');
    try
        idx = (size(X,2)-30*pSR):size(X,2); %select last 30 seconds of the data for plotting
        idx = idx(idx>0);
        if (pSR>2)
            pwelch(X(ChanFreq,idx),[],[],[],pSR)
        else
            pwelch(X(ChanFreq,idx))
        end
        s = sprintf('PSD of Channel %d Before Inteference Cancellation',ChanFreq);
        title(s);
    end
    
    
function MakeExample
    	fs = 500;
		n = 60*fs; %1-min sequence	
        m = 10; %number of channels
		t = 2*pi*(1:n)/fs;
		fline = 60 + randn; %ramdom interference frequency
		s = filter(1,[1,-0.99],100*randn(n,m))'; %1/f PSD
		p = bsxfun(@times,sin(fline*t+randn), (80+10*randn(m,1))) ...
         + bsxfun(@times,sin(2*fline*t+randn), (50+10*rand(m,1))) ...
		  + bsxfun(@times,sin(3*fline*t+randn), (20+10*randn(m,1))); % interference	
		gui_example = s + p;
        assignin('base','gui_example',gui_example);

        
        
    


% --- Executes just before removePLI_gui is made visible.
function removePLI_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to removePLI_gui (see VARARGIN)
MakeExample
% Choose default command line output for removePLI_gui
handles.output = hObject;


% TEXT annotations need an axes as parent so create an invisible axes which
% is as big as the figure
% Find all static text UICONTROLS whose 'Tag' starts with latex_
lbls = findobj(hObject,'-regexp','tag','latex_*');
allpt = [];
for i=1:length(lbls)
      l = lbls(i);
      % Get current text, position and tag
      pt = get(l,'parent');
      if i==1 
          handles.laxis(i) = axes('parent',pt,'units','normalized','position',[0 0 1 1],'visible','off');
          allpt = pt;
      else
          if any(allpt==pt)
            axes(handles.laxis(allpt==pt));
          else
              handles.laxis(i) = axes('parent',pt,'units','normalized','position',[0 0 1 1],'visible','off');
              allpt = [allpt pt];
          end
      end
      set(l,'units','normalized');
      s = get(l,'string');
      p = get(l,'position');
      t = get(l,'tag');
      % Remove the UICONTROL
      delete(l);
      % Replace it with a TEXT object 
      handles.(t) = text(p(1),p(2),s,'interpreter','latex');
      %set(handles.(t),'parent',pt);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes removePLI_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
removePLI_gui_InitFields(handles);



% --- Outputs from this function are returned to the command line.
function varargout = removePLI_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edtSRate_Callback(hObject, eventdata, handles)
% hObject    handle to edtSRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtSRate as text
%        str2double(get(hObject,'String')) returns contents of edtSRate as a double
SR = str2double(get(handles.edtSRate,'String'));
ACFreq = get(handles.popACFreq,'Value');
if ACFreq==2, 
    Harms = 1:floor(SR/2/50);
else
    Harms = 1:floor(SR/2/60);
end
if isempty(Harms) || any(isnan(Harms));
    set(handles.popHarm,'String','---');
else
    set(handles.popHarm,'String',Harms);
end



% --- Executes during object creation, after setting all properties.
function edtSRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtSRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popACFreq.
function popACFreq_Callback(hObject, eventdata, handles)
% hObject    handle to popACFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popACFreq contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popACFreq
SR = str2double(get(handles.edtSRate,'String'));
ACFreq = get(handles.popACFreq,'Value');
if ACFreq==2, 
    Harms = 1:floor(SR/2/50);
else
    Harms = 1:floor(SR/2/60);
end
if isempty(Harms) || any(isnan(Harms));
    set(handles.popHarm,'String','---');
else
    set(handles.popHarm,'String',Harms);
end


% --- Executes during object creation, after setting all properties.
function popACFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popACFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on popACFreq and none of its controls.
function popACFreq_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to popACFreq (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtB0_Callback(hObject, eventdata, handles)
% hObject    handle to edtB0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtB0 as text
%        str2double(get(hObject,'String')) returns contents of edtB0 as a double


% --- Executes during object creation, after setting all properties.
function edtB0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtB0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edtBinf_Callback(hObject, eventdata, handles)
% hObject    handle to edtBinf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtBinf as text
%        str2double(get(hObject,'String')) returns contents of edtBinf as a double


% --- Executes during object creation, after setting all properties.
function edtBinf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtBinf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtBst_Callback(hObject, eventdata, handles)
% hObject    handle to edtBst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtBst as text
%        str2double(get(hObject,'String')) returns contents of edtBst as a double


% --- Executes during object creation, after setting all properties.
function edtBst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtBst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtP0_Callback(hObject, eventdata, handles)
% hObject    handle to edtP0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtP0 as text
%        str2double(get(hObject,'String')) returns contents of edtP0 as a double


% --- Executes during object creation, after setting all properties.
function edtP0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtP0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtPinf_Callback(hObject, eventdata, handles)
% hObject    handle to edtPinf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPinf as text
%        str2double(get(hObject,'String')) returns contents of edtPinf as a double


% --- Executes during object creation, after setting all properties.
function edtPinf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPinf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtPst_Callback(hObject, eventdata, handles)
% hObject    handle to edtPst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPst as text
%        str2double(get(hObject,'String')) returns contents of edtPst as a double


% --- Executes during object creation, after setting all properties.
function edtPst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtW_Callback(hObject, eventdata, handles)
% hObject    handle to edtW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtW as text
%        str2double(get(hObject,'String')) returns contents of edtW as a double


% --- Executes during object creation, after setting all properties.
function edtW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popInputVar.
function popInputVar_Callback(hObject, eventdata, handles)
% hObject    handle to popInputVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

removePLI_gui_InitFields(handles);




% --- Executes during object creation, after setting all properties.
function popInputVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popInputVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtOutputVar_Callback(hObject, eventdata, handles)
% hObject    handle to edtOutputVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtOutputVar as text
%        str2double(get(hObject,'String')) returns contents of edtOutputVar as a double


% --- Executes during object creation, after setting all properties.
function edtOutputVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtOutputVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popHarm.
function popHarm_Callback(hObject, eventdata, handles)
% hObject    handle to popHarm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popHarm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popHarm


% --- Executes during object creation, after setting all properties.
function popHarm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popHarm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popChanFreq.
function popChanFreq_Callback(hObject, eventdata, handles)
% hObject    handle to popChanFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popChanFreq contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popChanFreq
InitPSD(handles)


% --- Executes during object creation, after setting all properties.
function popChanFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popChanFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popPlotChan.
function popPlotChan_Callback(hObject, eventdata, handles)
% hObject    handle to popPlotChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popPlotChan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popPlotChan
DoPlots(handles)

% --- Executes during object creation, after setting all properties.
function popPlotChan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popPlotChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnRemovePLI.
function btnRemovePLI_Callback(hObject, eventdata, handles)
% hObject    handle to btnRemovePLI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Run the Interference cancellation algorithm
try
    runRemovePLI(handles)
catch ME
    s = sprintf('Error occured during interference cancellation!\n%s\n Please check all the input fields and try again.',ME.message);
    msgbox(s ,'Error!','error','modal');
end
    


% --- Executes on button press in btnPaper.
function btnPaper_Callback(hObject, eventdata, handles)
% hObject    handle to btnPaper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
url = 'http://iopscience.iop.org/1741-2552/11/2/026017';
web(url,'-browser')

% --- Executes on button press in btnWebsite.
function btnWebsite_Callback(hObject, eventdata, handles)
% hObject    handle to btnWebsite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
url = 'https://github.com/mrezak/removePLI';
web(url,'-browser')


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text22.
function text22_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
url = 'https://github.com/mrezak/removePLI';
web(url,'-browser')



function edtSkip_Callback(hObject, eventdata, handles)
% hObject    handle to edtSkip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtSkip as text
%        str2double(get(hObject,'String')) returns contents of edtSkip as a double


% --- Executes during object creation, after setting all properties.
function edtSkip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtSkip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnPlot.
function btnPlot_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DoPlots(handles)
