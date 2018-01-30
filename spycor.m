function spycor

hObject = figure('Visible','off','Position',[260,300,800,600],...
    'toolbar','figure','MenuBar','none','NumberTitle','off','name','spycor');

handles = guidata(hObject);

% Choose default command line output for spycor
handles.output = hObject;

handles.current_directory = pwd; % nastaveni aktualniho adresare
handles.filepath = handles.current_directory;
handles.filename = '';

handles.N_spectra_limit = 3;  % The least number of spectra that can be treated.

handles.spectra_orig = zeros(2,2);
handles.spectra_corr = handles.spectra_orig;
handles.avg_spc = handles.spectra_orig;
handles.std_spc = zeros(size(handles.spectra_orig, 2), 1);
handles.corrIdx = cell(size(handles.spectra_orig, 2), 1);
handles.N_spectra = 1;
handles.recalculate = true(1, 1);
handles.calculated_once = zeros(1, 1);

handles.std_times = 30;
handles.std_times_separate = ones(1, 1) * handles.std_times;
handles.extra_del_points = 3;
handles.extra_del_points_separate = ones(1, 1) * handles.extra_del_points;
handles.extent = 0;
handles.extent_separate = cell(1, 1);
handles.extent_separate(:) = {handles.extent};
handles.range = [1, 1];  % indices of selected range
handles.extent_separate = cell(1, 1);
handles.extent_separate(:) = {handles.range};
handles.all = true;
handles.all_used = false;  % if all spectra had the same values of settings and has not been set to idividual settings
% yet, the default value for them will be the same as the settings for all spectra.
handles.chosen_spectrum = 1;

% functional handles
h=spycor_functions;
[handles.spycor_load, handles.find_corrIdx, handles.vec2str]=h{:};

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

handles.range_pushbutton = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [.55 .15 .115 .8],...
    'String', 'range', 'Parent', handles.settings_panel,...
    'TooltipString', 'Sets used spectral range',...
    'FontSize', 10, 'Callback', @range_pushbutton_Callback);
handles.extent_pushbutton = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [.68 .15 .115 .8],...
    'String', 'extent', 'Parent', handles.settings_panel,...
    'TooltipString', 'Sets which neighbour spectra are used for spike detection',...
    'FontSize', 10, 'Callback', @extent_pushbutton_Callback);

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
  if size(data_orig, 2) < handles.N_spectra_limit + 1
    h_errordlg=errordlg(sprintf('The file must contain %d spectra at least!', handles.N_spectra_limit));
    waitfor(h_errordlg);
    return;
  end
  handles.x_scale = data_orig(:,1);
  handles.spectra_orig = data_orig(:,2:end);
  N = size(handles.spectra_orig, 2);
  handles.N_spectra = N;

  handles.spectra_corr = handles.spectra_orig;
  handles.avg_spc = handles.spectra_orig;
  handles.std_spc = zeros(N, 1);
  handles.corrIdx = cell(N, 1);

  handles.chosen_spectrum = 1;
  handles.recalculate = true(N, 1);
  handles.calculated_once = zeros(N, 1);

  handles.std_times_separate = ones(N, 1) * handles.std_times;
  handles.extra_del_points_separate = ones(N, 1) * handles.extra_del_points;
  handles.extent = 0;
  handles.extent_separate = cell(N, 1);
  handles.extent_separate(:) = {handles.extent};
  handles.range = [1, length(handles.x_scale)];
  handles.range_separate = cell(N, 1);
  handles.range_separate(:) = {handles.range};

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
if handles.chosen_spectrum < handles.N_spectra
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
    handles.recalculate = true(handles.N_spectra, 1);
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
    handles.recalculate = true(handles.N_spectra, 1);
else
    handles.extra_del_points_separate(ii) = str2double(get(hObject, 'String'));
    handles.recalculate(ii) = true;
end

guidata(hObject, handles);

treat_data(hObject, ii);

plot_function(hObject);


% --- Executes on button press in range_pushbutton.
function range_pushbutton_Callback(hObject, eventdata, defans)
% hObject    handle to range_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% defams     defaut answer
handles = guidata(hObject);
ii = handles.chosen_spectrum;
N = handles.N_spectra;

if nargin < 3
    if handles.all
        defans_var = handles.range;
    else
        defans_var = handles.range_separate{ii};
    end
    if isempty(defans_var)
        defans = {''};
    else
        defans = sprintf('%d, %d', handles.x_scale(defans_var(1)), handles.x_scale(defans_var(2)));
        for jj = 3:2:length(defans_var)
            defans = sprintf('%s, %d, %d', defans, handles.x_scale(defans_var(jj)),...
                handles.x_scale(defans_var(jj + 1)));
        end
        defans = {defans};
    end
end

prompt = {sprintf(['Specify x range of treatment. Beware that this drops the results from\n'...
                   'previous range! The range should be specified as the pairs of limits\n'...
                   'of desired subranges, for example 500, 700, 1230, 1600 means\n'...
                   '<500, 700> U <1230, 1600>. Empty field excludes spectrum from treatment.\n'...
                   'The spectral limits are: <%d, %d>'],...
    handles.x_scale(1), handles.x_scale(end))};
dlg_title = 'Selection of spectral range';
num_lines = 1;
options.Resize = 'on';
range_text = inputdlg(prompt, dlg_title, num_lines, defans, options);

% Test cancel
if isequal(range_text,{})
    return;
end

% Test input validity
try
    range_vals = eval(sprintf('[%s]', range_text{1}));
    if ~isempty(range_vals) && (~isnumeric(range_vals) || ~isvector(range_vals) || rem(length(range_vals), 2))
        msgID = 'range:BadInput';
        msg = 'Range must be even number of numbers.';
        exception = MException(msgID, msg);
        throw(exception)
    end
    % convert ranges to indices and test if they contain at least one point
    range_inds = [];
    for jj = 1:2:length(range_vals)
        x = handles.x_scale(handles.x_scale > range_vals(jj) & handles.x_scale < range_vals(jj + 1));
        if ~isempty(x)
            [Y, I1] = min(abs(handles.x_scale - range_vals(jj)));
            [Y, I2] = min(abs(handles.x_scale - range_vals(jj + 1)));
            range_inds = [range_inds, I1:I2];
        end
    end
    range_inds = unique(range_inds, 'sorted');
    if (~isempty(range_vals) && isempty(range_inds))
        msgID = 'range:BadInput';
        exception = MException(msgID, 'In this range aren''t any spectral points. Please, specify correct range.');
        throw(exception)
    end
catch ME
    getReport(ME)
    h_errordlg=errordlg(sprintf('%s', ME.message), 'Input error');
    waitfor(h_errordlg);
    range_pushbutton_Callback(hObject, eventdata, range_text)
    return
end

if ~isempty(range_inds)
    range = range_inds(1);
    rangeinds = 1:(length(range_inds) - 1);
    rangeinds = rangeinds(diff(range_inds)~=1);
    rangeinner = [range_inds(rangeinds); range_inds(rangeinds + 1)];
    range = [range, rangeinner(:)', range_inds(end)];
else
    range = [];
end

% save range
if handles.all
    if isequal(range, handles.range)
        return;
    end
    handles.recalculate = true(N, 1);
    handles.range = range;
else
    if isequal(range, handles.range_separate{ii})
        return;
    end
    handles.recalculate(ii) = true;
    handles.range_separate{ii} = range;
end

guidata(hObject, handles);
change_spectrum(hObject);


% --- Executes on button press in extent_pushbutton.
function extent_pushbutton_Callback(hObject, eventdata, defans)
% hObject    handle to extent_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% defams     defaut answer
handles = guidata(hObject);

ii = handles.chosen_spectrum;
N_limit = handles.N_spectra_limit;
N = handles.N_spectra;

if nargin < 3
    if handles.all
        defans_var = handles.extent;
    else
        defans_var = handles.extent_separate{ii};
    end
    if length(defans_var) == 1 && defans_var == 0
        defans = {num2str(defans_var)};
    else
        % choose default answer to shorter one from range with zero or without it
        defans1 = handles.vec2str(defans_var);
        defans2 = handles.vec2str(sort([defans_var, 0]));
        if length(defans1) > length(defans2)
            defans = {defans2};
        else
            defans = {defans1};
        end
    end
end

prompt = {sprintf(['Specify range of surrounding spectra, which is used to detect spikes.\n'...
                   'Use MatLab vector creation notation, for example, if you want to use\n'...
                   '3 spectra from the left and 2 spectra from the right, use "-3:2", for\n'...
                   'spectra -5, -3, -2, 1, 2, 3 use "[-5, -3, -2, 1:3]". If you want to use\n'...
                   'all the other spectra, use "0". Edge cases will be cropped, so be sure\n'...
                   'that there are at least %d surrounding spectra used for each spectrum.'],...
    N_limit - 1)};
dlg_title = 'Selection of surrounding spectra';
num_lines = 1;
options.Resize = 'on';
extent_text = inputdlg(prompt, dlg_title, num_lines, defans, options);

% Test cancel
if isequal(extent_text,{})
    return;
end

% Test input validity
try
    extent = eval(extent_text{1});
    if isempty(extent) || ~isnumeric(extent) || ~isvector(extent) || any(rem(extent, 1))
        msgID = 'extent:BadInput';
        msg = 'Extent must be MatLab vector of integers.';
        exception = MException(msgID, msg);
        throw(exception)
    end
    if length(extent) > 1 || extent ~= 0
        extent = sort(setdiff(extent, 0));  % make extent sorted set of unique values
        % reduce extent to only valid surrounding spectra and test if minimal number of spectra is greater then N_limit
        if handles.all
            n_min = 2 * (N - 1);  % initialize to maximal number of surrounding spectra on both sides
            extent = extent((extent > -N) & (extent < N));  % reduce extent only to relevant numbers
            if length(extent) == 2 * (N - 1)  % test if extend contains all possibly surrounding spectra
                extent = 0;
            else  % find minimal overlap of extend with spectral set
                for jj = 1:N
                    extent_length = length(extent((extent + jj > 0) & (extent + jj <= N)));
                    if extent_length < n_min
                        n_min = extent_length;
                    end
                end
            end
            msg = sprintf(['Extent must have at least %d surrounding spectra for all treated spectra '...
                '(including corner cases).'], N_limit - 1);
        else  % not handles.all
            extent = extent((extent + ii > 0) & (extent + ii <= N));  % reduce extent only to relevant numbers
            n_min = length(extent);
            if n_min == N - 1  % test if extend contains all possibly surrounding spectra
                extent = 0;
            end
            msg = sprintf('Extent must have at least %d surrounding spectra.', N_limit - 1);
        end
        if (n_min < N_limit - 1)
            msgID = 'extent:BadInput';
            exception = MException(msgID, msg);
            throw(exception)
        end
    end
catch ME
    getReport(ME)
    h_errordlg=errordlg(sprintf('%s', ME.message), 'Input error');
    waitfor(h_errordlg);
    extent_pushbutton_Callback(hObject, eventdata, extent_text)
    return
end

% save extent
if handles.all
    if isequal(extent, handles.extent)
        return;
    end
    handles.recalculate = true(N, 1);
    handles.extent = extent;
else
    if isequal(extent, handles.extent_separate{ii})
        return;
    end
    handles.recalculate(ii) = true;
    handles.extent_separate{ii} = extent;
end

guidata(hObject, handles);
change_spectrum(hObject);


% --- Executes on button press in all_checkbox.
function all_checkbox_Callback(hObject, eventdata)
% hObject    handle to all_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
handles = guidata(hObject);

handles.all = get(hObject, 'Value');

N = handles.N_spectra;
ii = handles.chosen_spectrum;

if handles.all
    handles.recalculate =...
        ~handles.calculated_once...
        | (handles.std_times_separate ~= handles.std_times * ones(N, 1))...
        | (handles.extra_del_points_separate ~= handles.extra_del_points * ones(N, 1));
    for jj = 1:N
        handles.recalculate(jj) = handles.recalculate(jj) | ~isequal(handles.extent_separate{jj}, handles.extent)...
             | ~isequal(handles.range_separate{jj}, handles.range);
    end

    set(handles.all_pushbutton, 'Enable', 'off')
else
    if ~handles.all_used
        handles.std_times_separate = ones(N, 1) * handles.std_times;
        handles.extra_del_points_separate = ones(N, 1) * handles.extra_del_points;
        handles.extent_separate(:) = {handles.extent};
        handles.range_separate(:) = {handles.range};
    end
    handles.recalculate(ii) =...
        (handles.std_times_separate(ii) ~= handles.std_times)...
        | (handles.extra_del_points_separate(ii) ~= handles.extra_del_points)...
        | ~isequal(handles.extent_separate{ii}, handles.extent)...
        | ~isequal(handles.range_separate{ii}, handles.range);

    set(handles.all_pushbutton, 'Enable', 'on')
end

handles.all_used = true;

guidata(hObject, handles);

change_spectrum(hObject);


% --- Executes on button press in all_pushbutton.
function all_pushbutton_Callback(hObject, eventdata)
% hObject    handle to all_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
handles = guidata(hObject);

ii = handles.chosen_spectrum;

handles.recalculate =...
    (handles.std_times_separate ~= handles.std_times_separate(ii))...
    | (handles.extra_del_points_separate ~= handles.extra_del_points_separate(ii));
for jj = 1:handles.N_spectra
    if ~isequal(handles.extent_separate{jj}, handles.extent_separate{ii})
        handles.recalculate(jj) = true;
        handles.extent_separate{jj} = handles.extent_separate{ii};
    end
    if ~isequal(handles.range_separate{jj}, handles.range_separate{ii})
        handles.recalculate(jj) = true;
        handles.range_separate{jj} = handles.range_separate{ii};
    end
end
handles.std_times_separate = ones(handles.N_spectra, 1) * handles.std_times_separate(ii);
handles.extra_del_points_separate = ones(handles.N_spectra, 1) * handles.extra_del_points_separate(ii);
guidata(hObject, handles);


function treat_data(hObject, spc_no)
handles = guidata(hObject);

ii = spc_no;

handles.calculated_once(ii) = true;

N = handles.N_spectra;
handles.spectra_corr(:,ii) = handles.spectra_orig(:,ii);

if handles.all
    std_times = handles.std_times;
    extra_del_points = handles.extra_del_points;
    if handles.extent
        extent = handles.extent + ii;
    else
        extent = [1:(ii - 1), (ii + 1):N];
    end
    range = handles.range;
else
    std_times = handles.std_times_separate(ii);
    extra_del_points = handles.extra_del_points_separate(ii);
    if handles.extent
        extent = handles.extent_separate{ii} + ii;
    else
        extent = [1:(ii - 1), (ii + 1):N];
    end
    range = handles.range_separate{ii};
end
if isempty(range)
    handles.corrIdx{ii} = {};
    handles.std_spc(ii) = 0;
    handles.avg_spc(:,ii) = handles.spectra_orig(:,ii);
else
    inds = [];
    for jj = 1:2:length(range)
        inds = [inds, range(jj):range(jj + 1)];
    end
    dx = (handles.x_scale(2) - handles.x_scale(1)) * 1.5;
    [corrIdx, handles.avg_spc(inds,ii), handles.std_spc(ii)] = ...
        handles.find_corrIdx(handles.x_scale(inds), handles.spectra_orig(inds,:), dx, ii, std_times, extent);
    if ~isempty(corrIdx) % controll, if there is any spike in
        % subtracting identified spikes from the spectrum
        for jj = 1:size(corrIdx, 2)
            corrIdx{jj} = inds(corrIdx{jj});  % convert corrIdx from inds indexes to general indexes
            I1 = corrIdx{jj}(1);
            range_no = find(range(1:2:end) < I1, 1, 'last');  % identify range number
            if I1 > range(2 * range_no - 1) + extra_del_points - 1
                idx1 = I1 - extra_del_points;
            else
                idx1 = range(2 * range_no - 1);
            end
            I2 = corrIdx{jj}(end);  % convert corrIdx from inds indexes to general indexes
            if I2 < range(2 * range_no) - extra_del_points + 1
                idx2 = I2 + extra_del_points;
            else
                idx2 = range(2 * range_no);
            end
            handles.spectra_corr(idx1:idx2,ii) = handles.avg_spc(idx1:idx2,ii);
        end
    end
    handles.corrIdx{ii} = corrIdx;
end
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
    range = handles.range;
else
    std_times = handles.std_times_separate(ii);
    range = handles.range_separate{ii};
end
if isempty(range)
    plot(handles.x_scale, handles.avg_spc(:,ii), '-', 'Color', [192, 192, 192] / 255);
else
    if range(1) ~= 1
        cinds = 1:range(1);
        plot(handles.x_scale(cinds), handles.avg_spc(cinds, ii), '-', 'Color', [192, 192, 192] / 255);
        hold on;
    end
    for jj = 3:2:length(range)
        cinds = range(jj - 1):range(jj);
        plot(handles.x_scale(cinds), handles.avg_spc(cinds, ii), '-', 'Color', [192, 192, 192] / 255);
        hold on;
    end
    if range(end) ~= length(handles.x_scale)
        cinds = range(end):length(handles.x_scale);
        plot(handles.x_scale(cinds), handles.avg_spc(cinds, ii), '-', 'Color', [192, 192, 192] / 255);
        hold on;
    end

    inds = range(1):range(2);
    pH(1) = plot(handles.x_scale(inds), handles.avg_spc(inds,ii),'-k');
    hold on;
    plot(handles.x_scale(inds), handles.avg_spc(inds,ii) + std_times * handles.std_spc(ii), ':k');
    pH(2) = plot(handles.x_scale(inds), handles.spectra_orig(inds,ii),'-b');
    for jj = 3:2:length(range)
        inds = range(jj):range(jj + 1);
        plot(handles.x_scale(inds), handles.avg_spc(inds,ii),'-k');
        plot(handles.x_scale(inds), handles.avg_spc(inds,ii) + std_times * handles.std_spc(ii), ':k');
        plot(handles.x_scale(inds), handles.spectra_orig(inds,ii),'-b');
    end

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
            pH(3) = plot(handles.x_scale(idx1:idx2), handles.spectra_orig(idx1:idx2,ii), '-r', 'LineWidth', 2);
            pH(4) = plot(handles.x_scale(idx1:idx2), handles.spectra_corr(idx1:idx2,ii), '-g', 'LineWidth', 2);
        end
        legend(pH, {'spc', 'avg', 'spike', 'corrected'}, 'Location', 'best');
    else
        legend(pH, {'spc', 'avg'}, 'Location', 'best');
    end
end

xmin = min(handles.x_scale);
xmax = max(handles.x_scale);
dx = (xmax - xmin) / 50;
axis([xmin - dx, xmax + dx, -inf, inf])

guidata(hObject,handles);
