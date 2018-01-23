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
handles.avg_spc = handles.spectra_orig;
handles.std_spc = zeros(size(handles.spectra_orig, 2), 1);
handles.corrIdx = cell(size(handles.spectra_orig, 2), 1);
handles.N_spectra = 1;
handles.recalculate = ones(1, 1);
handles.calculated_once = zeros(1, 1);

handles.std_times = 30;
handles.std_times_separate = ones(1, 1) * handles.std_times;
handles.extra_del_points = 3;
handles.extra_del_points_separate = ones(1, 1) * handles.extra_del_points;
handles.all = true;
handles.all_used = false;  % if all spectra had the same values of settings and has not been set to idividual settings
% yet, the default value for them will be the same as the settings for all spectra.
handles.chosen_spectrum = 1;

% functional handles
h=spycor_functions;
[handles.spycor_load, handles.find_corrIdx]=h{:};

handles.file_menu = uimenu('Label','File');
handles.load_menuitem = uimenu(handles.file_menu, 'Label','Load','Accelerator','O',...
    'Callback',@load_menuitem_Callback);
handles.load_menuitem = uimenu(handles.file_menu, 'Label','Save','Accelerator','S',...
    'Callback',@save_menuitem_Callback);

handles.select_spec_panel = uipanel('Title','Select spectrum','Units',...
    'normalized','Position',[.05,.88,.35,.1],'Visible','off');
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
    'normalized','Position',[.42,.88,.53,.1],'Visible','off');
handles.times_text = uicontrol('Style','text','Units','normalized','Position',...
    [.01 .35 .13 .4],'String','times:','Parent',handles.settings_panel,...
    'TooltipString','Multiples of standard deviation to resolve spike','FontSize',10);
handles.times_edit = uicontrol('Style','edit','Units','normalized','Position',...
    [.14 .15 .10 .8],'String',num2str(handles.std_times),'Parent',handles.settings_panel,...
    'TooltipString','Multiples of standard deviation to resolve spike','FontSize',10,...
    'Callback', @times_edit_Callback, 'CreateFcn', @times_edit_CreateFcn);

handles.extra_del_points_text = uicontrol('Style','text','Units','normalized','Position',...
    [.25 .35 .12 .4],'String','extra:','Parent',handles.settings_panel,...
    'TooltipString','Number of points to be deleted around spike on both sides.','FontSize',10);
handles.extra_del_points_edit = uicontrol('Style','edit','Units','normalized','Position',...
    [.37 .15 .10 .8],'String',num2str(handles.extra_del_points),'Parent',handles.settings_panel,...
    'TooltipString','Number of points to be deleted around spike on both sides.','FontSize',10,...
    'Callback', @extra_del_points_edit_Callback, 'CreateFcn', @extra_del_points_edit_CreateFcn);

if handles.all
    all_pushbutton_enable = 'off';
else
    all_pushbutton_enable = 'on';
end
handles.all_pushbutton = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [.81 .15 .085 .8],...
    'String', 'all', 'Parent', handles.settings_panel,...
    'TooltipString', 'Sets the current settings for all spectra',...
    'Enable', all_pushbutton_enable, 'FontSize', 10, 'Callback', @all_pushbutton_Callback);
handles.all_checkbox = uicontrol('Style', 'checkbox', 'Units', 'normalized', 'Position', [.90 .35 .10 .4],...
    'String', 'all', 'Parent', handles.settings_panel,...
    'Value', handles.all,...
    'Tooltipstring', ['All spectra have the same settings. If set to true, the settings for all spectra will be set'...
    ' to settings for the current spectrum. Custom settings will be restored if it is set to false again.'],...
    'Callback', @all_checkbox_Callback);

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
  N = size(handles.spectra_orig, 2);
  handles.N_spectra = N;

  handles.spectra_corr = handles.spectra_orig;
  handles.avg_spc = handles.spectra_orig;
  handles.std_spc = zeros(N, 1);
  handles.corrIdx = cell(N, 1);

  handles.chosen_spectrum = 1;
  handles.recalculate = ones(N, 1);
  handles.calculated_once = zeros(N, 1);

  handles.std_times_separate = ones(handles.N_spectra, 1) * handles.std_times;
  handles.extra_del_points_separate = ones(handles.N_spectra, 1) * handles.extra_del_points;


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
  
  treat_data(hObject, handles.chosen_spectrum);
  
  change_spectrum(hObject);
end

function save_menuitem_Callback(hObject, eventdata)
% hObject    handle to save_menuitem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

if any(handles.recalculate)
    I = find(handles.recalculate);
    for ii = 1:length(I)
        treat_data(hObject, I(ii));
    end
end

handles = guidata(hObject);

savedata = [handles.x_scale ,handles.spectra_corr];
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
        handles.chosen_spectrum = num_chos_spec;

        guidata(hObject, handles);
        change_spectrum(hObject);
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

guidata(hObject, handles);
change_spectrum(hObject);


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

guidata(hObject, handles);
change_spectrum(hObject);


% --- Executes during object creation, after setting all properties.
function times_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to times_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function times_edit_Callback(hObject, eventdata)
% hObject    handle to times_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of times_edit as text
%        str2double(get(hObject,'String')) returns contents of times_edit as a double

handles = guidata(hObject);

ii = handles.chosen_spectrum;

if handles.all
    handles.std_times = str2double(get(hObject, 'String'));
    handles.recalculate = ones(handles.N_spectra, 1);
else
    handles.std_times_separate(ii) = str2double(get(hObject, 'String'));
    handles.recalculate(ii) = true;
end

guidata(hObject, handles);
  
treat_data(hObject, ii);

plot_function(hObject);

% --- Executes during object creation, after setting all properties.
function extra_del_points_edit_CreateFcn(hObject, eventdata)
% hObject    handle to extra_del_points_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function extra_del_points_edit_Callback(hObject, eventdata)
% hObject    handle to extra_del_points_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

handles = guidata(hObject);

ii = handles.chosen_spectrum;

if handles.all
    handles.extra_del_points = str2double(get(hObject, 'String'));
    handles.recalculate = ones(handles.N_spectra, 1);
else
    handles.extra_del_points_separate(ii) = str2double(get(hObject, 'String'));
    handles.recalculate(ii) = true;
end

guidata(hObject, handles);

treat_data(hObject, ii);

plot_function(hObject);


% --- Executes on button press in all_checkbox.
function all_checkbox_Callback(hObject, eventdata)
% hObject    handle to all_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
handles = guidata(hObject);

handles.all = get(hObject, 'Value');

N = handles.N_spectra;
ii = handles.chosen_spectrum;

if ~handles.all_used && ~handles.all
    handles.std_times_separate = ones(N, 1) * handles.std_times;
    handles.extra_del_points_separate = ones(N, 1) * handles.extra_del_points;
end

if handles.all
    handles.recalculate =...
        ~handles.calculated_once...
        | (handles.std_times_separate ~= handles.std_times * ones(N, 1))...
        | (handles.extra_del_points_separate ~= handles.extra_del_points * ones(N, 1));

    set(handles.all_pushbutton, 'Enable', 'off')
else
    handles.recalculate(ii) =...
        (handles.std_times_separate(ii) ~= handles.std_times)...
        | (handles.extra_del_points_separate(ii) ~= handles.extra_del_points);

    set(handles.all_pushbutton, 'Enable', 'on')
end

handles.all_used = true;

guidata(hObject, handles);

change_spectrum(hObject);


% --- Executes on button press in choose_spec_pushbutton.
function all_pushbutton_Callback(hObject, eventdata)
% hObject    handle to all_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
handles = guidata(hObject);

ii = handles.chosen_spectrum;

handles.recalculate =...
    (handles.std_times_separate ~= handles.std_times_separate(ii))...
    | (handles.extra_del_points_separate ~= handles.extra_del_points_separate(ii));
handles.std_times_separate = ones(handles.N_spectra, 1) * handles.std_times_separate(ii);
handles.extra_del_points_separate = ones(handles.N_spectra, 1) * handles.extra_del_points_separate(ii);
guidata(hObject, handles);


function treat_data(hObject, spc_no)

handles = guidata(hObject);

ii = spc_no;

handles.calculated_once(ii) = true;

N = handles.N_spectra;

if handles.all
    std_times = handles.std_times;
    extra_del_points = handles.extra_del_points;
else
    std_times = handles.std_times_separate(ii);
    extra_del_points = handles.extra_del_points_separate(ii);
end

[corrIdx, handles.avg_spc(:,ii), handles.std_spc(ii)] = ...
        handles.find_corrIdx(handles.x_scale, handles.spectra_orig, ii, std_times);

if ~isempty(corrIdx) % controll, if there is any spike in
    % subtracting identified spikes from the spectrum
    handles.spectra_corr(:,ii) = handles.spectra_orig(:,ii);
    for jj = 1:size(corrIdx, 2)
        if corrIdx{jj}(1) > extra_del_points
            idx1 = corrIdx{jj}(1) - extra_del_points;
        else
            idx1 = 1;
        end
        if corrIdx{jj}(end) < length(handles.x_scale) - extra_del_points + 1
            idx2 = corrIdx{jj}(end) + extra_del_points;
        else
            idx2 = length(handles.x_scale);
        end
        handles.spectra_corr(idx1:idx2,ii) = handles.avg_spc(idx1:idx2,ii);
    end
end

handles.corrIdx{ii} = corrIdx;
handles.recalculate(ii) = false;

guidata(hObject,handles);


function change_spectrum(hObject)

handles = guidata(hObject);

ii = handles.chosen_spectrum;

set(handles.chosen_spec_text, 'String', sprintf('Chosen spectrum: %d', ii), ...
    'BackgroundColor', 'green');
if handles.all
    set(handles.times_edit, 'String', num2str(handles.std_times));
    set(handles.extra_del_points_edit, 'String', num2str(handles.extra_del_points));
else
    set(handles.times_edit, 'String', num2str(handles.std_times_separate(ii)));
    set(handles.extra_del_points_edit, 'String', num2str(handles.extra_del_points_separate(ii)));
end

guidata(hObject, handles);

if handles.recalculate(ii)
    treat_data(hObject, ii);
end

plot_function(hObject);



function plot_function(hObject)

handles = guidata(hObject);

ii = handles.chosen_spectrum;

cla(handles.main_axes);

if handles.all
    std_times = handles.std_times;
else
    std_times = handles.std_times_separate(ii);
end

pH(1) = plot(handles.x_scale, handles.avg_spc(:,ii),'-k');
hold on;
plot(handles.x_scale, handles.avg_spc(:,ii) + std_times * handles.std_spc(ii),':k');

pH(end + 1) = plot(handles.x_scale, handles.spectra_orig(:,ii),'-b');
if ~isempty(handles.corrIdx{ii})
    for jj = 1:size(handles.corrIdx{ii},2)
        if handles.corrIdx{ii}{jj}(1) > handles.extra_del_points
            idx1 = handles.corrIdx{ii}{jj}(1)-handles.extra_del_points;
        else
            idx1 = 1;
        end
        if handles.corrIdx{ii}{jj}(end) < length(handles.x_scale) - handles.extra_del_points + 1
            idx2 = handles.corrIdx{ii}{jj}(end)+handles.extra_del_points;
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

xmin = min(handles.x_scale);
xmax = max(handles.x_scale);
dx = (xmax - xmin) / 50;
axis([xmin - dx, xmax + dx, -inf, inf])

guidata(hObject,handles);
