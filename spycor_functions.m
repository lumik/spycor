function h=textrapol_functions

h={@spycor_load, @find_corrIdx};

function [status,filename,file_path,spectra]=spycor_load(do_test,start_directory,...
FilterSpec,dialog_title,multiselect)
%--------------------------------------------------------------------------
% Nacteni spekter ze souboru (spektra ve sloupcich, 1. sloupec je x-ova skala) 
%--------------------------------------------------------------------------
% Syntaxe funkce:
% [status,filename,file_path,spectra]=readdata(do_test,start_directory,...
% FilterSpec,dialog_title,multiselect)
%--------------------------------------------------------------------------
% Vstupni parametry:
% do_test -> promenna do_test nabyva hodnoty 0,1. Pro hodnotu 1 se testuje,
% jestli je x-ova skala setridena (viz poznamka).
% start_directory -> pocatecni adresar, ktery se vypise v dialogovem okne
% FilterSpec -> 'Cell Array' pro specifikaci typu souboru pro otevreni (tak 
% jak je definovano ve funkci Matlabu 'uigetfile') 
% dialog_title -> retezec, ktery se vypisuje v dialogovem okne 
% multiselect -> 1: vybirame vice souboru, 0: lze vybrat pouze jeden soubor
%--------------------------------------------------------------------------
% Vystupni parametry:
% status -> v pripade stisku cancel pri nacitani dat ma status hodnotu 0.
% V opacnem pripade ma status hodnotu 1.
% filename -> jmeno nacteneho souboru (v pripade, ze se vybira vice souboru
% v rezimu multiselect=1, jsou jmena nactenych souboru v promenne typu Cell
% array).
% file_path -> absolutni cesta k souboru
% spectra -> sloupce v matici spectra odpovidaji jednotlivym nactenym spektrum 
% (prvni sloupec v matici je x-ova skala). v pripade stisku cancel pri
% nacitani dat je spectra=[]
%--------------------------------------------------------------------------
% Poznamka: 1) Pokud prvni sloupec dat (x-ove hodnoty) neni setriden (coz muze
% znacit chybu), pak je nabidnuto jeho vzestupne setrideni. Podle tohoho
% prvniho sloupce se pak setridi ostatni sloupce (spektralni intenzity). S
% takto setridenymi spektry se lepe pracuje pri vykreslovani grafu,
% orezavani, apod.
% 2)  Nacitana Spektra mohou byt ulozena v textovem souboru nebo v binarnim
% souboru (format matlabu: pripona mat). 
%--------------------------------------------------------------------------

% Vyber souboru:
[filename,file_path] = uigetfile(FilterSpec,dialog_title,...
fullfile(start_directory,'*.*'),'Multiselect',multiselect);
if (isequal(filename,0) || isequal(file_path,0)) % stisknuto cancel
 status=0;
 spectra=[];
else
 status=1;   
 spec=strcat(file_path,filename); % specifikace pro nasledne nacteni dat ze souboru
 try
  spectra=load(spec); % nacteni spekter ze souboru (textovy nebo binarni soubor "mat").
 catch
  status=0;
  filename=0;
  file_path=0;
  spectra=[];
  return
 end
 extension=filename((strfind(filename,'.')+1):1:end); % pripona souboru
 if strcmp(extension,'mat') % Spektra se nactou z "mat" souboru 
  polozky_spektra=fieldnames(spectra);   
  spectra=eval(['spectra.' polozky_spektra{1}]); % v souboru "mat" je ulozena
  % pouze jedna promenna (pozadovana spektra)
 end
 if do_test==1 &&  ~issorted(spectra(:,1)) && ~issorted(flipdim(spectra(:,1),1))
  % x-ove hodnoty nejsou setrideny ani vzestupne, ani sestupne -> muze 
  % signalizovat chybu (je nabidnuta moznost spektra vzestupne setridit)
  vypis=['Spectra have unsorted x-values (probably due to the error). '...
  'Do you want to sort them in increasing order ?'];    
  tridit = questdlg(vypis,'Loading data','Yes','No','Yes');
  if strcmp(tridit,'Yes')
   [spectra(:,1),indexy]=sort(spectra(:,1)); % tridi se x-ovy hodnoty
   for i=2:1:size(spectra,2) % y-ovy hodnoty pro kazdy spektrum se tridi podle x
    y=spectra(:,i);
    spectra(:,i)=y(indexy);
   end        
  end
 end
end

function [corrIdx, avg_y, std_y] = find_corrIdx(x, y, times)
avg_y = y;
N = size(y, 2);
std_y = zeros(N,1);
corrIdx = cell(N,1);

for ii = 1:N
    idx = setdiff(1:N,ii);
    avg_y(:,ii) = mean(y(:,idx), 2);
    std_y(ii) = mean(std(y(:,idx), 0, 2));

    corrIdxVec = find(y(:, ii) > avg_y(:,ii) + times * std_y(ii));
    kk = 1;
    dx = (x(2) - x(1)) * 1.5;
    corrIdx{ii} = {};
    if ~isempty(corrIdxVec)
        corrIdx{ii}{kk} = [];
        for jj = 1:length(corrIdxVec)
            corrIdx{ii}{kk} = [corrIdx{ii}{kk}; corrIdxVec(jj)];
            if (jj < length(corrIdxVec)) && (x(corrIdxVec(jj+1)) - x(corrIdxVec(jj)) > dx)
                kk = kk + 1;
                corrIdx{ii}{kk} = [];
            end
        end
    end
end