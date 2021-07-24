function [status,result]=systemnomatlab(command)
    % Trims Matlab libraries from the PATH variable on a system call
    %
    % MATLAB has the habit of prepending itself to the path ruining all kinds of
    % dependencies during system calls with very hard to trace back error
    % messages.
    % This function is called as like the matlab "system" but the PATH variable
    % is edited to not have any path with the matlab prepended.
    %
    % May cause a failure if you Need some of those matlab libraries later in
    % the path, especially if Matlab is actually first on your path, then
    % probably you shouldn't need this anyway.

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
