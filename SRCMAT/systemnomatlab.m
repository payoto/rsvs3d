function [status,result]=systemnomatlab(command)
    
    [~, matlabedPath]=system('echo %PATH%');
    
    if isunix
        pathSep = ':';
    else
        pathSep = ';';
    end
    
    pathElms = regexp(matlabedPath,pathSep,'split');
    if isunix
        pathTrimFunc = @(str) ['PATH//:', str, ':/:'];
        pathTrimStart = 'PATH=:$PATH: && ';
        pathTrimEnd = ' && PATH=${PATH#:}; PATH=${PATH%:} && ';
    else
        pathTrimFunc = @(str) [str, pathSep];
        pathTrimStart = 'set path=%PATH:';
        pathTrimEnd = '=% & ';
    end
    
    kk=1;
    pathTrimCmd = pathTrimStart;
    while kk<=numel(pathElms) && ~isempty(regexpi(pathElms{kk},'matlab'))
        pathTrimCmd=[pathTrimCmd, pathTrimFunc(pathElms{kk})];
        kk=kk+1;
    end
    pathTrimCmd=[pathTrimCmd, pathTrimEnd];
    
    [status, result] = system([pathTrimCmd, command]);

end