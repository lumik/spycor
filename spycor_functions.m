function h=spycor_functions

h={@spycor_load, @find_corrIdx, @vec2str};

function [status, filename, file_path, spectra] = spycor_load(do_test, start_directory, filter_spec, dialog_title,...
    multiselect)
% Load spectra from file (spectra are in columns, 1st column is x axis)
% [status, filename, file_path, spectra] =...
%     spycor_load(do_test, start_directory, filter_spec, dialog_title,...
%     multiselect)
% inputs:
% do_test          true test if the x scale is sorted.
% start_directory  starting directory of dialog window.
% filter_spec      'Cell Array' of filetype specifications as defined in
%                  'uigetfile' MatLab function.
% dialog_title     title of load file dialog.
% multiselect      if true more files is possible to select.
% outputs:
% status           0 if cancel was pushed, 1 otherwise.
% filename         name of selected file (if multiselect == 1, Cell array of
%                  filenames).
% file_path        absolutni path to file
% spectra          spectra in column matrix, 1st column is x scale. Empty
%                  matrix if cancel was pushed.
% Notes: 1) If the first column (x scale) is not sorted (it can be an error),
%           sorting is offered. For the treatment spectra must be sorted.
%        2) Loaded spectra can be in matlab binary file (must have mat
%           extension) or text file (all the other extensions).

% file selection
[filename, file_path] = uigetfile(filter_spec, dialog_title, fullfile(start_directory, '*.*'),...
    'Multiselect', multiselect);

status = 0;
if (isequal(filename, 0) || isequal(file_path, 0))  % pushed cancel
    spectra = [];
    return
end

full_filename = strcat(file_path, filename);  % selected full file name

try
    spectra=load(full_filename);  % load spectra (autodetects extension).
catch
    status = 0;
    filename = 0;
    file_path = 0;
    spectra = [];
    return
end

extension = filename((find(filename == '.', 1, 'last') + 1):1:end);  % file extension

if strcmp(extension, 'mat')  % binary data from '.mat' file
    data_fields = fieldnames(spectra);
    spectra = eval(['spectra.' data_fields{1}]); % In 'mat' file can be only one variable (matrix with spectra)
end

if do_test == 1 &&  ~issorted(spectra(:,1))
    % x scale is not sorted
    msg = ['Spectra have unsorted x-values (probably due to the error). '...
        'Do you want to sort them in ascending order?'];
    sort_answer = questdlg(msg, 'Loading data', 'Yes', 'No', 'Yes');
    if strcmp(sort_answer, 'Yes')
        [spectra(:,1), ind] = sort(spectra(:,1));
        spectra(:,2:end) = spectra(ind,2:end);
    else
        status = 0;
        filename = 0;
        file_path = 0;
        spectra = [];
        return
    end
end
status = 1;


function [corrIdx, avg_y, std_y] = find_corrIdx(x, y, dx, ii, times, extent)
% [corrIdx, avg_y, std_y] = find_corrIdx(x, y, dx, ii, times, extent)
% Finds groups of adjacent indices, which are more deviated from the average
% than the times * mean(stdev in the corresponding values of the surrounding spectra (y) specified by extent).
N = size(y, 2);

% find indices of spectra for average calculation
idx = setdiff(extent(extent > 0 & extent <= N), ii);

avg_y = mean(y(:,idx), 2);
std_y = mean(std(y(:,idx), 0, 2));

corrIdxVec = find(y(:,ii) > avg_y + times * std_y);
kk = 1;
corrIdx = {};
if ~isempty(corrIdxVec)
    corrIdx{kk} = [];
    for jj = 1:length(corrIdxVec)
        corrIdx{kk} = [corrIdx{kk}; corrIdxVec(jj)];
        if (jj < length(corrIdxVec)) && (x(corrIdxVec(jj + 1)) - x(corrIdxVec(jj)) > dx)
            kk = kk + 1;
            corrIdx{kk} = [];
        end
    end
end


function strvec = vec2str(vec)
% Converts vector to string

NO = 3;  % minimal number of elements for x:y notation
DNO = 4;  % minimal number of elements for x:y:z notation
N = length(vec);
if N == 1
    strvec = num2str(vec);
    return;
else
    brackets = false;
    ino = 0;  % counter for x:y notation
    idno = 0;  % counter for x:y:z notation
    d = vec(2) - vec(1);  % difference
    strvec = num2str(vec(1));
    for ii = 2:N
        % if the difference from previous is one
        if vec(ii) - vec(ii - 1) == 1 && idno == 0
            ino = ino + 1;
        % if the difference from previous is d
        elseif vec(ii) - vec(ii - 1) == d
            idno = idno + 1;
        % leap
        else
            brackets = true;
            % if continous sequence satisfactory long
            if ino >= NO
                strvec = sprintf('%s:%d', strvec, vec(ii - 1));
            % if continous sequence with difference satisfactory long
            elseif idno >= DNO
                strvec = sprintf('%s:%d:%d', strvec, d, vec(ii - 1));
            % if continous sequence is not satisfactory long
            else
                % continous sequence
                if ino > 0
                    count = ino;
                % other cases
                else
                    count = idno;
                end
                for jj = (ii - count):ii - 1
                    strvec = sprintf('%s, %d', strvec, vec(jj));
                end
            end
            strvec = sprintf('%s, %d', strvec, vec(ii));
            if ii < N
                d = vec(ii + 1) - vec(ii);
            end
            ino = 0;
            idno = 0;
        end
    end
    % treat the last iteration
    % if continous sequence satisfactory long
    if ino >= NO
        strvec = sprintf('%s:%d', strvec, vec(N));
    % if continous sequence with difference satisfactory long
    elseif idno >= DNO
        strvec = sprintf('%s:%d:%d', strvec, d, vec(N));
    % if continous sequence is not satisfactory long
    else
        brackets = true;
        % continous sequence
        if ino > 0
            count = ino;
        % other cases
        else
            count = idno;
        end
        for jj = (N - count):N
            strvec = sprintf('%s, %d', strvec, vec(jj));
        end
    end
end

if brackets
    strvec = sprintf('[%s]', strvec);
end
