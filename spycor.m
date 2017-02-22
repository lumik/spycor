function spycor

hObject = figure('Visible','off','Position',[260,300,800,600],...
    'toolbar','figure','MenuBar','none','NumberTitle','off','name','spycor');

handles = guidata(hObject);

% Choose default command line output for spycor
handles.output = hObject;

handles.current_directory = pwd; % nastaveni aktualniho adresare
handles.filepath = handles.current_directory;
handles.filename = '';

handles.spectra_orig = zeros(2,2);
handles.spectra_corr = handles.spectra_orig;
handles.N_spectra = 1;

handles.std_times = 30;
handles.around_spike = 3;
handles.chosen_spectrum = 1;

handles.minN = 2; % number of spectra, where disappears spike when spike
% in the current spectrum is subtracted

% functional handles
h=spycor_functions;
[handles.spycor_load, handles.find_corrIdx]=h{:};

handles.file_menu = uimenu('Label','File');
handles.load_menuitem = uimenu(handles.file_menu, 'Label','Load','Accelerator','O',...
    'Callback',@load_menuitem_Callback);
handles.load_menuitem = uimenu(handles.file_menu, 'Label','Save','Accelerator','S',...
    'Callback',@save_menuitem_Callback);

handles.select_spec_panel = uipanel('Title','Select spectrum','Units',...
    'normalized','Position',[.05,.88,.4,.1],'Visible','off');
handles.chosen_spec_text = uicontrol('Style','pushbutton','Units','normalized',...
    'Position',[.47 .15 .5 .8],'String','Chosen spectrum:','Parent',...
    handles.select_spec_panel,'BackgroundColor','green','Enable',...
    'inactive');
handles.choose_spec_pushbutton = uicontrol('Style', 'pushbutton','Units','normalized','Position',...
    [.02 .15 .13 .8],'String','#','Parent',handles.select_spec_panel,...
    'TooltipString','Will chose one spectrum','Visible','on',...
    'Enable','on','FontSize',10, 'Callback', @choose_spec_pushbutton_Callback);
handles.down_pushbutton = uicontrol('Style','pushbutton','Units','normalized','Position',...
    [.17 .15 .13 .8],'String','<','Parent',handles.select_spec_panel,...
    'TooltipString','Will chose one spectrum','FontSize',10,...
    'Callback', @down_pushbutton_Callback);
handles.up_pushbutton = uicontrol('Style','pushbutton','Units','normalized','Position',...
    [.32 .15 .13 .8],'String','>','Parent',handles.select_spec_panel,...
    'TooltipString','Will chose one spectrum','FontSize',10,...
    'Callback', @up_pushbutton_Callback);
handles.main_axes = axes('Units','normalized','Position',[0.05,0.05,0.9,0.8]);

handles.settings_panel = uipanel('Title','Settings','Units',...
    'normalized','Position',[.55,.88,.4,.1],'Visible','off');
handles.times_text = uicontrol('Style','text','Units','normalized','Position',...
    [.02 .15 .17 .8],'String','times:','Parent',handles.settings_panel,...
    'TooltipString','Multiples of standard deviation to resolve spike','FontSize',10);
handles.times_edit = uicontrol('Style','edit','Units','normalized','Position',...
    [.21 .15 .13 .8],'String',num2str(handles.std_times),'Parent',handles.settings_panel,...
    'TooltipString','Multiples of standard deviation to resolve spike','FontSize',10,...
    'Callback', @times_edit_Callback, 'CreateFcn', @times_edit_CreateFcn);

handles.around_text = uicontrol('Style','text','Units','normalized','Position',...
    [.36 .15 .17 .8],'String','around:','Parent',handles.settings_panel,...
    'TooltipString','Number of points to be deleted around spike','FontSize',10);
handles.around_edit = uicontrol('Style','edit','Units','normalized','Position',...
    [.55 .15 .13 .8],'String',num2str(handles.around_spike),'Parent',handles.settings_panel,...
    'TooltipString','Number of points to be deleted spike','FontSize',10,...
    'Callback', @around_edit_Callback, 'CreateFcn', @around_edit_CreateFcn);

%Make the UI visible.
set(hObject,'Visible','on');

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function load_menuitem_Callback(hObject, eventdata)
% hObject    handle to load_menuitem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = guidata(hObject);
handles = guidata(hObject);
filter_spec={'*.*','All Files (*.*)';
 '*.mat','MAT-files (*.mat)';
 '*.txt;*.dat','text files (*.txt,*.dat)'}; 
[status,soubor,cesta,data_orig]=handles.spycor_load(0,handles.current_directory,...
  filter_spec,'Load Data','Off');
if status==0 % Pokud dojde k chybe pri nacitani dat (stisknuti cancel),
  % zobrazi se chybova hlaska.   
  h_errordlg=errordlg('No data load or invalid format of data !',...
      'loading data');
  waitfor(h_errordlg);
else % V pripade nacteni dat se provadi vse dalsi v Callbacku (prepisou se
  % predchozi data a smazou predchozi fity).

  handles.x_scale = data_orig(:,1);
  handles.spectra_orig = data_orig(:,2:end);
  handles.spectra_corr = handles.spectra_orig;
  handles.N_spectra = size(handles.spectra_orig,2);
  handles.chosen_spectrum = 1;

  set(handles.select_spec_panel, 'Visible','on');
  set(handles.settings_panel, 'Visible','on');
  set(handles.chosen_spec_text, 'String', ...
      sprintf('Chosen spectrum: %d', handles.chosen_spectrum), ...
      'BackgroundColor', 'green');
  set(handles.down_pushbutton,'Enable','off');
  if(handles.N_spectra == 1)
      set(handles.up_pushbutton,'Enable','off');
  else
      set(handles.up_pushbutton,'Enable','on');
  end
  handles.filepath = cesta;
  handles.current_directory = cesta;
  handles.filename = soubor;
  
  guidata(hObject, handles);
  
  treat_data(hObject);
  
  plot_function(hObject);
end

function save_menuitem_Callback(hObject, eventdata)
% hObject    handle to save_menuitem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

savedata = [handles.x_scale,handles.spectra_corr];
I = find(handles.filename == '.', 1, 'last');
filename = [handles.filepath,handles.filename(1:I-1),'_kor.txt'];
fprintf('Saving file:\n%s\n',filename);
save(filename,'savedata','-ascii');
fprintf('Done.\n');

% --- Executes on button press in choose_spec_pushbutton.
function choose_spec_pushbutton_Callback(hObject, eventdata)
% hObject    handle to choose_spec_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
prompt = {sprintf('Select one spectrum\n(integer from 1 to %d)',...
    handles.N_spectra)}; 
dlg_title_matrix = 'Selection of spectrum';
num_lines = 1;
options.Resize='on';
odpoved_num_of_spec=inputdlg(prompt,dlg_title_matrix,num_lines,{'',''},...
    options);
%--------------------------------------------------------------------------
% Testovani, zda nebylo stisknuto cancel
%--------------------------------------------------------------------------
if ~isequal(odpoved_num_of_spec,{})
    [num_chos_spec,status]=str2num(odpoved_num_of_spec{1});
    if num_chos_spec < 1 || num_chos_spec > handles.N_spectra || ...
            status == 0 || ...
            ~isequal(floor(num_chos_spec), num_chos_spec) || ...
            size(num_chos_spec,2) ~= 1
        % Nepripustne hodnoty
        set(handles.chosen_spec_text,'String', 'Incorrect number', ...
        'BackgroundColor','red'); % Chybna hodnota a jeji cervena indikace
        %errordlg('You selected incorret number of spectra','Incorret input','on'); 
    else % Pripustne hodnoty mezi matice se ulozi do globalnich promennych
        if num_chos_spec > 1
            set(handles.down_pushbutton, 'Enable', 'on');
        else
           set(handles.down_pushbutton, 'Enable', 'off');
        end
        if num_chos_spec<handles.N_spectra
            set(handles.up_pushbutton, 'Enable', 'on');
        else
            set(handles.up_pushbutton, 'Enable', 'off');
        end
        set(handles.chosen_spec_text, 'String', ...
            sprintf('Chosen spectrum: %d', num_chos_spec), ...
            'BackgroundColor', 'green');
        handles.chosen_spectrum = num_chos_spec;
        guidata(hObject, handles);
        plot_function(hObject);
    end
end 

% --- Executes on button press in down_pushbutton.
function down_pushbutton_Callback(hObject, eventdata)
% hObject    handle to down_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if handles.chosen_spectrum <= 2
    handles.chosen_spectrum = 1;
    set(handles.down_pushbutton,'Enable','off');
else
    handles.chosen_spectrum = handles.chosen_spectrum - 1;
end
if handles.chosen_spectrum < handles.N_spectra;
    set(handles.up_pushbutton,'Enable','on');
end

set(handles.chosen_spec_text, 'String', ...
    sprintf('Chosen spectrum: %d', handles.chosen_spectrum), ...
    'BackgroundColor', 'green');
guidata(hObject,handles);
plot_function(hObject);

% --- Executes on button press in up_pushbutton.
function up_pushbutton_Callback(hObject, eventdata)
% hObject    handle to up_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if handles.chosen_spectrum >= handles.N_spectra - 1
    handles.chosen_spectrum = handles.N_spectra;
    set(handles.up_pushbutton,'Enable','off');
else
    handles.chosen_spectrum = handles.chosen_spectrum + 1;
end
if handles.chosen_spectrum > 1
    set(handles.down_pushbutton,'Enable','on');
end

set(handles.chosen_spec_text, 'String', ...
    sprintf('Chosen spectrum: %d', handles.chosen_spectrum), ...
    'BackgroundColor', 'green');
guidata(hObject,handles);
plot_function(hObject);

% --- Executes during object creation, after setting all properties.
function times_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kak_slider_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function times_edit_Callback(hObject, eventdata)
% hObject    handle to kak_slider_step_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of times_edit as text
%        str2double(get(hObject,'String')) returns contents of times_edit as a double

handles = guidata(hObject);

handles.std_times = str2double(get(hObject, 'String'));

guidata(hObject, handles);
  
treat_data(hObject);

plot_function(hObject);

% --- Executes during object creation, after setting all properties.
function around_edit_CreateFcn(hObject, eventdata)
% hObject    handle to kak_slider_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function around_edit_Callback(hObject, eventdata)
% hObject    handle to kak_slider_step_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of times_edit as text
%        str2double(get(hObject,'String')) returns contents of times_edit as a double

handles = guidata(hObject);

handles.around_spike = str2double(get(hObject, 'String'));

guidata(hObject, handles);
  
treat_data(hObject);

plot_function(hObject);

function treat_data(hObject)

handles = guidata(hObject);

[corrIdx, handles.avg_spc, handles.std_spc] = ...
    handles.find_corrIdx(handles.x_scale, handles.spectra_orig, handles.std_times);

% correction of false spikes caused by influence of spike on the average
% spectrum

% iteration over the spectra
for ii = 1:handles.N_spectra
    if ~isempty(corrIdx{ii}) % controll, if there is any spike in
        % subtracting identified spikes from the spectrum
        corrSpc = handles.spectra_orig;
        for jj = 1:size(corrIdx{ii},2)
            if corrIdx{ii}{jj}(1) > handles.around_spike
                idx1 = corrIdx{ii}{jj}(1)-handles.around_spike;
            else
                idx1 = 1;
            end
            if corrIdx{ii}{jj}(end) < length(handles.x_scale) - handles.around_spike + 1
                idx2 = corrIdx{ii}{jj}(end)+handles.around_spike;
            else
                idx2 = length(handles.x_scale);
            end
            corrSpc(idx1:idx2,ii) = handles.avg_spc(idx1:idx2,ii);
        end
        
        % finding spikes in the data set with ii-th spectrum corrected for
        % spikes
        corrIdx2 = ...
            handles.find_corrIdx(handles.x_scale, corrSpc, handles.std_times);
        
        otherIdx = setdiff(1:handles.N_spectra,ii);
        
        % iteration on all spikes
        for jj = 1:size(corrIdx{ii},2)
            currN = 0; % number of spectra, which contained spike before
            % correction of spectrum ii but do not contain it after
            % correction of this spectrum
            kkll = [];
            
            % iteration over all the other spectra
            for kk = otherIdx
                if ~isempty(corrIdx{kk}) % control, if there is any spike
                    % in this spectrum
                    success = false; % true, if there is a spike in the
                    % spectrum kk on the same position as the spike jj in
                    % the spectrum ii
                    
                    ll = 1; % iterates oves all spikes in the other spectrum
                    kkll2 = [];
%                     while ll <= size(corrIdx{kk},2) && ~success
                    while ll <= size(corrIdx{kk},2)
                        if ~isempty(intersect(corrIdx{kk}{ll}, corrIdx{ii}{jj}))
                            success = true;
                            kkll2 = [kkll2; kk, ll]; % index of occurence of spike
                            % at the same position in the spectrum kk
                            % as in the spectrum ii
                        end
                        ll = ll + 1;
                    end
                    
                    % if the spike jj of the spectrum ii was found also
                    % in the spectrum kk
                    if success
                        % and is not found in the spectrum kk after
                        % subtraction of spikes from the spectrum ii
                        % increase the number of false spikes currN
                        % and save the index of the spike in the kk
                        % spectrum
                        if isempty(corrIdx2{kk})
                            kkll = [kkll; kkll2];
                            currN = currN + 1;
                        else
                            ll = 1;
                            falsespike = true;
                            
                            while ll <= size(corrIdx2{kk},2) && falsespike
                                if ~isempty(intersect(corrIdx{kk}{ll}, corrIdx{ii}{jj}))
                                    falsespike = false;
                                end
                                ll = ll + 1;
                            end
                            if falsespike
                                kkll = [kkll; kkll2];
                                currN = currN + 1;
                            end
                        end
                    end
                end
            end
            % if there is detected false spike in more than or equal to
            % the threshold handles.minN
            if currN >= handles.minN
                % remove them from the other spectra
                mm = 1;
                while mm <= size(kkll, 1)
                    notSpikeIdx = 1:size(corrIdx{kkll(mm,1)},2);
                    currkk = kkll(mm,1);
                    samespike = true;
                    while samespike
                        notSpikeIdx = setdiff(notSpikeIdx,kkll(mm,2));
                        if mm < size(kkll, 1) && currkk == kkll(mm+1,1)
                            mm = mm + 1;
                        else
                            samespike = false;
                        end
                    end
                    tempCorrIdx = corrIdx{kkll(mm,1)};
                    corrIdx{kkll(mm,1)} = {};
                    if ~isempty(notSpikeIdx)
                        for nn = 1:length(notSpikeIdx)
                            corrIdx{kkll(mm,1)}{nn} = tempCorrIdx{notSpikeIdx(nn)};
                        end
                    end
                    mm = mm + 1;
                end
            end
        end
    end
end

% construction of corrected spectra
handles.spectra_corr = handles.spectra_orig;

for ii = 1:handles.N_spectra
    if ~isempty(corrIdx{ii}) % controll, if there is any spike in
        for jj = 1:size(corrIdx{ii},2)
            if corrIdx{ii}{jj}(1) > handles.around_spike
                idx1 = corrIdx{ii}{jj}(1)-handles.around_spike;
            else
                idx1 = 1;
            end
            if corrIdx{ii}{jj}(end) < length(handles.x_scale) - handles.around_spike + 1
                idx2 = corrIdx{ii}{jj}(end)+handles.around_spike;
            else
                idx2 = length(handles.x_scale);
            end
            handles.spectra_corr(idx1:idx2,ii) = handles.avg_spc(idx1:idx2,ii);
        end
    end
end

handles.corrIdx = corrIdx;

guidata(hObject,handles);

function plot_function(hObject)

handles = guidata(hObject);

cla(handles.main_axes);

ii = handles.chosen_spectrum;
pH(1) = plot(handles.x_scale, handles.avg_spc(:,ii),'-k');
hold on;
plot(handles.x_scale, handles.avg_spc(:,ii)+handles.std_times*handles.std_spc(ii),':k');

pH(end + 1) = plot(handles.x_scale, handles.spectra_orig(:,ii),'-b');
if ~isempty(handles.corrIdx{ii})
    for jj = 1:size(handles.corrIdx{ii},2)
        if handles.corrIdx{ii}{jj}(1) > handles.around_spike
            idx1 = handles.corrIdx{ii}{jj}(1)-handles.around_spike;
        else
            idx1 = 1;
        end
        if handles.corrIdx{ii}{jj}(end) < length(handles.x_scale) - handles.around_spike + 1
            idx2 = handles.corrIdx{ii}{jj}(end)+handles.around_spike;
        else
            idx2 = length(handles.x_scale);
        end
        pH(3) = plot(handles.x_scale(idx1:idx2), handles.spectra_orig(idx1:idx2,ii),'-r','LineWidth',2);
        pH(4) = plot(handles.x_scale(idx1:idx2), handles.spectra_corr(idx1:idx2,ii),'-g','LineWidth',2);
    end
    legend(pH, {'spc','avg','spike','corrected'},'Location','best');
else
    legend(pH, {'spc','avg'},'Location','best');
end

aLims = axis;
xmin = min(handles.x_scale);
xmax = max(handles.x_scale);
dx = (xmax - xmin) / 50;
axis([xmin - dx, xmax + dx, aLims(3:4)])

guidata(hObject,handles);
