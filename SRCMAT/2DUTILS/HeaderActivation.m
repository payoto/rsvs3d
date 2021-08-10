%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2016
%
%          Function for Creation of Header
%                     Folders
%                Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [] = HeaderActivation(funcHandles,funcDir)
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.

    c=evalc('display(funcHandles)');
    pattern='@\w*';
    fprintf('\nWRITING AUTOMATIC FUNCTIONS')
    funcHandlesNamesCell=regexp(c,pattern,'match');
    funcDir=MakePathCompliant(funcDir);
    if ~isdir(funcDir)

        mkdir(funcDir);
    end
    jj=0;
    callerFunc=dbstack(1);
    for ii=1:length(funcHandlesNamesCell)
        funcName=funcHandlesNamesCell{ii}(2:end);
        funcHandleVarname=[funcHandlesNamesCell{ii}(2:end),'_Handle'];
        eval(['global ',funcHandleVarname])
        eval([funcHandleVarname,'=funcHandles{ii};']);
        if  ~isFindFuncFile(funcDir,[funcName,'.m'])
            WriteContainerFunctionFile(funcName,funcDir,callerFunc(1).name);
        else
            jj=jj+1;
        end
    end
    addpath(funcDir);
    if jj~=0
        fprintf(' - %i of %i files skipped\n',jj,ii);
    end
end

function pathName=MakePathCompliant(pathName)

    compStr=computer;
    if strcmp(compStr(1:2),'PC')
        pathName=regexprep(pathName,'/','\\');
    else

        pathName=regexprep(pathName,'\\','/');
    end
end

function []=WriteContainerFunctionFile(funcName,funcDir,callerName)
    fID=-1;
    t=now;
    Dt=0;
    while(fID<0 && Dt<0.01)
        fID=fopen(MakePathCompliant([funcDir,'\',funcName,'.m']),'w');
        Dt=now-t;
    end
    funcText{1}=['function [varargout]=',funcName,'(varargin)'];
    funcText{end+1}=['%% ',callerName];
    funcText{end+1}=['global ',funcName,'_Handle'];
    funcText{end+1}=['nOut=nargout(',funcName,'_Handle',');'];
    funcText{end+1}=['nOutReq=nargout;'];
    funcText{end+1}=['nOut(nOut<0)=nOutReq;'];
    funcText{end+1}=['[varargout{1:nOut}]=',funcName,'_Handle','(varargin{:});'];
    funcText{end+1}=['end'];

    for ii=1:length(funcText)

        fprintf(fID,[funcText{ii},'\n']);

    end
    fclose('all');

end

function [isFound]=isFindFuncFile(rootDir,strDir)

    subDir=dir(rootDir);
    subDir(1:2)=[];
    ii=1;
    isFound=false;
    while ~isFound && ii<=numel(subDir)

        isFound=strcmp(subDir(ii).name,strDir);
        ii=ii+1;

    end


end
