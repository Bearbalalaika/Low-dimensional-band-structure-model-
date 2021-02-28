%import selected data files from directory
%MB 18/05/2006
%ver 1.01
%when exits 
%leaves following variables apart from imported data:
%dataFileNames, dataPathName, dataN 

if exist('dataFileNames') || exist('dataPathName') || exist('dataN')
    answer=questdlg('Clear old variables: dataFileNames, dataPathName, dataN, dataLog ?', ...
        'Clear variables ?','Yes','Cancel','Yes' );
    if answer=='Yes'
        %get list of vars and delete them
        list=evalin('base','dataFileNames');
        if  iscell(list)
            h=waitbar(0,'Clearing...');  
            set( get( get(h,'Children'),'Title') ,'Interpreter','none');
            n=size(list,2);
            for i=1:n
                str=char(list(1,i));
                %function str=guessVarName(str, n) from ple.m
                %if first char is a letter its ok - matlab did not modified it
                %if it is a number, matlab added a X to a variable name
                str=fliplr(char(str));
                str=fliplr(str(1,strfind(str,'.')+1:end));
                if length(str2num(str(1,1)))
                    str=sprintf('clear X%s',str);%-------------------------
                else
                    str=sprintf('clear %s',str);%-------------------------
                end
                evalin('base',str);
                waitbar(i/n,h,str);
            end
            close(h);
        end
        %delete list
        eval('clear dataFileNames dataPathName dataN dataLog');
    else
        break;
    end;
end

str=eval('pwd');
if strcmp(str,'/Applications/MATLAB71')
    cd '/Users/marcin/Documents/data';
end

%get dir listing with mask
[dataFileNames, dataPathName, dataFilterIndex]=uigetfile({'*.xy';'*.dat';'*.txt';'*.*'},'Step 1/2: Select Data Files','MultiSelect','on');
if isequal(dataFileNames,0)
    break;
end
cd (dataPathName);

[logName, logPathName]=uigetfile({'*.log';'*.txt';'*.*'},'Step 2/2: Select Log File','MultiSelect','off');
if isequal(logName,0)
    %
else
    dataLog='';
    %load text log file into dataLog, line by line
    fid = fopen(logName, 'rt');
    y = 0;
    while feof(fid) == 0
        y=y+1;
        tline = fgetl(fid);
        dataLog{y}=tline;
    end
    fclose(fid);
end

%load data
dataN=size(dataFileNames,2);
h=waitbar(0,'Importing...'); 
set( get( get(h,'Children'),'Title') ,'Interpreter','none');

for n=1:dataN; 
    name=char(dataFileNames(1,n));
    load(name,'-ascii');
    str=strcat('Importing...',name, ' done');
    waitbar(n/dataN,h,str);
end
close(h);


clear ans i n h list name answer str dataFilterIndex;
clear fid y tline logName logPathName;

[file, path]=uiputfile('*.mat');
if file(1,1)~=0
    path=strcat(path,file);
    clear file;
    save(path);
end
clear path;
clear file;

