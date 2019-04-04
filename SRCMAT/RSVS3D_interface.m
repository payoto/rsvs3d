function [] = RSVS3D_interface()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
    include_Utilities
end


%% Parameter reading and parsing

% Read default parameters
function [defaultparam]=RSVS3D_DefaultParam(pathRSVS, defaultName)
    
    if~exist('pathRSVS', 'var'); pathRSVS = '.';end
    if~exist('defaultName', 'var')
        defaultName = './defaultconfigRSVS3D.json';
    end
    
    if isunix()
        fileExt = '';
    else
        fileExt = '.exe';
    end
    
    pathDefaultConf = [pathRSVS, filesep, 'RSVS3D', fileExt];
    command =[pathDefaultConf, ' --noexec=', defaultName];
    [s,o] = systemnomatlab(command);
    if s
        errorstruct.identifier='rsvs3D:param:default';
        errorstruct.message=['Failure to process command "',command,'"'];
        error(errorstruct)
    end
    [defaultparam]=RSVS3D_ReadParam(defaultName);
end

function [paramOut]=RSVS3D_ReadParam(fileIn)
    
    fid = fopen(fileIn, 'r');
    string = '';
    while ~feof(fid)
        string = [string, squeeze(fgetl(fid))];
    end
    fclose(fid);
    paramOut = jsondecode(string);
    [paramOut.structdat]=GetStructureData(paramOut);
    
end

% Write new parameter files
function [success]=RSVS3D_WriteParam(paramIn, fileOut)
    
    if(numel(ExtractVariables({'inputpoints'},paramIn))<2)
        paramIn.grid.voronoi=rmfield(paramIn.grid.voronoi,'inputpoints');
    end
    
    fid = fopen(fileOut ,'w');
    strJson = jsonencode(rmfield(paramIn,'structdat'));
    success=fprintf(fid, '%s', strJson);
    fclose(fid);
    
end
% Handle a few standard modifications

%% Callers for the RSVS3D

% Caller to handle correcting the file structure
% This caller sets the outputs to a specific directory
function [param3D]=RSVS3D_CallInDirectory(param3D, targDir, varargin)
    
    if exist('targDir', 'var') && ~isempty(targDir)
        param3D.files.ioout.outdir = targDir;
    end
end

% Caller to generate meshes
% Achieved by calling RSVS3D with 0 steps and no exports
function [param3D]=RSVS3D_CallMeshes(param3D, varargin)
    
    param3D.snak.maxsteps = 0;
    [param3D]=RSVS3D_CallInDirectory(param3D, varargin{:});
end


% Caller to run export
% Run in full "loading" mode and 0 steps
function [param3D]=RSVS3D_CallExport(param3D, exportConf,voluMesh, snakeMesh, ...
        snakeFile, varargin)
    if(~isempty(exportConf))
        for ii = 1:numel(exportConf)
            if(numel(exportConf{ii})~=2)
                errstruct.message = sprintf(...
                    ["The number of 'exportconfig' parameters should be even",...
                    " (%i)"], numel(exportConf{ii}));
                errstruct.identifier = 'rsvs3D:caller:exportConfig:Odd';
                error(errstruct);
            end
        end
        param3D.files.exportconfig = exportConf;
    end
    param3D.snak.maxsteps = 0;
    if (isempty(voluMesh) || isempty(snakeFile) || isempty(snakeMesh))
        errstruct.message = sprintf(...
            ["Not all mesh files defined: voluMesh (%i), voluMesh (%i), ",...
            " snakeFile (%i)"], isempty(voluMesh), isempty(snakeFile),...
            isempty(snakeMesh));
        errstruct.identifier = 'rsvs3D:caller:exportConfig:NoMesh';
        error(errstruct);
    end
    [param3D]=RSVS3D_CallLoadMeshes(param3D, voluMesh, snakeMesh, ...
        snakeFile,varargin{:});
end

% caller to run RSVS process
% Full loading
function [param3D]=RSVS3D_CallLoadMeshes(param3D,voluMesh, snakeMesh, ...
        snakeFile,varargin)
    
    if ~exist('voluMesh', 'var'); voluMesh = ''; end
    if ~exist('snakeMesh', 'var'); snakeMesh = ''; end
    if ~exist('snakeFile', 'var'); snakeFile = ''; end
    
    if (isempty(snakeMesh) || isempty(voluMesh))
        errstruct.message = sprintf(['Cannot do a partial mesh load; \n', ...
            'Volume mesh: "%s"\nSnake mesh: "%s"'], voluMesh, snakeMesh);
        errstruct.identifier = 'rsvs3D:caller:partialMeshLoad';
        error(errstruct)
    end
    
    param3D.grid.activegrid = 'load';
    param3D.files.ioin.volumeshname = voluMesh;
    param3D.files.ioin.snakemeshname = snakeMesh;
    param3D.files.ioin.snakefile= snakeFile;
    [param3D]=RSVS3D_CallInDirectory(param3D, varargin{:});
end

% Base caller
function [outRSVS3D]=RSVS3D_Run(param3D, configOutput)
    % RSVS3D_CALLER Call the RSVS 3D executable with a predefined config.
    %  PARAM3D: The set of 3D-RSVS parameters to be called
    %  CONFIGOUTPUT: A structure configuring calling options for the RSVS3D
    %    The default structure can be built using RSVS3D_CALLCONFIG
    
    if ~exist('configOutput','var'); configOutput = RSVS3D_CallConfig();end
    
    % Write parameter file
    writeStatus=RSVS3D_WriteParam(param3D, configOutput.paramOut);
    
    if writeStatus<=0 % Not enough characters written
        errstruct.identifier = 'rsvs3D:caller:writeparam';
        errstruct.message = 'Parameter file was not correctly written';
        error(errstruct)
    end
    
    % Call RSVS
    [cmd]=RSVS3D_Config2Cmd(configOutput);
    [status, output]=systemnomatlab(cmd);
    
    if status~=0 % Program should end with status==0
        errstruct.identifier = 'rsvs3D:caller:execution';
        fmtMessage = ['Exit status of "RSVS3D" was not 0 (%i); \n', ...
            ' - command: \n%s\n - output:\n%s\n'];
        errstruct.message = sprintf(fmtMessage, status, cmd, output);
        error(errstruct)
    end
    
    % Check results for output folder
    [outRSVS3D, errorstatus] = RSVS3D_ProcessCall(output);
    if errorstatus
        warning('Errors detected while processic the RSVS3D output.')
    end
end

function [defaultConfig]=RSVS3D_CallConfig()
    % Sets where to call the executable from.
    
    % Path to the executable
    defaultConfig.executablePath = 'SRCC\RSVS3D';
    % Execution path, path from which this command should be executed
    defaultConfig.executionPath = '.';
    % Commands to prepend and post after the RSVS3D call
    defaultConfig.preCmd = '';
    defaultConfig.postCmd = '';
    % unique temp file
    if ~exist('temp','dir'), mkdir([cd,filesep,'temp']);end
    defaultConfig.paramOut = [tempname([cd,filesep,'temp']),'_3dparam.json'];
    defaultConfig.deleteParam = false;
    
end

function [cmd]=RSVS3D_Config2Cmd(config)
    % Builds the command from the config structure
    % structure generated by "RSVS3D_CallConfig"
    cmd = config.preCmd;
    if cmd; cmd = [cmd, ' && '];end
    cmd = [cmd, 'cd ', config.executionPath, ' && '];
    cmd = [cmd, config.executablePath, ' -l "', config.paramOut,'"'];
end

function [outRSVS3D, errorstatus] = RSVS3D_ProcessCall(charOut)
    % postprocess a call to RSVS3D
    
    outRSVS3D.char = regexp(charOut,'\n','split');
    errorstatus = 0;
    % Find errors
    errorStrings = {'error', 'unusual'};
    posError = RSVS3D_ProcessFindLines(outRSVS3D.char, errorStrings);
    errorstatus = numel(posError);
    % Find convergence
    convString = 'conv: ';
    posConv = RSVS3D_ProcessFindLines(outRSVS3D.char, convString);
    convstruct = repmat(RSVS3D_convstruct(),[1 0]);
    for i = posConv
        [convstruct(i)]=RSVS3D_ProcessConvergence(outRSVS3D.char{i},...
            convString);
    end
    outRSVS3D.conv = convstruct;
    % Find output directory
    dirString = 'Output folder:';
    posDir = RSVS3D_ProcessFindLines(outRSVS3D.char, dirString);
    if numel(posDir)>1
        warning('More than one folder output detected')
    end
    outRSVS3D.folder=RSVS3D_ProcessCallFolder(outRSVS3D.char{posDir(1)}, ...
        dirString);
end

function folder=RSVS3D_ProcessCallFolder(line, dirString)
    dirCells = regexp(squeeze(line), dirString, 'split');
    if (numel(dirCells)<1)
        warning('Output processing failed, returning multiple output folders.');
    elseif (numel(dirCells)==0)
        error('Output processing failed, no output folder found');
    end
    folder=dirCells{2};
end

function convstruct = RSVS3D_convstruct()
    convstruct.volume = 0;
    convstruct.velocity = 0;
end

function [convstruct]=RSVS3D_ProcessConvergence(line, convString)
    convstruct = RSVS3D_convstruct();
    
    fmtString = [convString, '(vol) %f (vel) %f'];
    posMatch = regexp(line, convString);
    [out, nOut, err] = sscanf(line(posMatch(1):end), fmtString);
    
    if nOut~=2
        error(err)
    end
    convstruct.volume = out(1);
    convstruct.velocity = out(2);
end

function [pos]=RSVS3D_ProcessFindLines(charAll, pattern)
    func = @(charAll, pat) find(~(cellfun(@isempty,regexp(charAll, pat))));
    D = size(charAll);
    pos = [];
    if iscell(pattern)
        for i = numel(pattern)
            posTemp = func(charAll, pattern{i});
            pos = [pos, posTemp];
        end
    else
        pos = func(charAll, pattern);
    end
    
end

%% callers for SU2

function  [su2meshparam] = RSVS3D_SU2Config()
    % Loads the configuration for SU2
    % Need loading of the tetgen::apiparam class for SU2 configuration
    
    % XYZ CFD domain lower bound
    su2meshparam.lowerB = -[15 15 15];
    % XYZ CFD domain upper bound
    su2meshparam.upperB = [15 15 15];
    % Edge lengths at regular intervals
    su2meshparam.edgelengths = [0.0250, 0.2000, 0.5000, 0, 0, 0, 2.0000];
    % Edge length at which two snake points are considered in the same 
    % place
    su2meshparam.distanceTol = 1e-3;
    % Surface edge lengths at minimum curvature and max curvature
    su2meshparam.surfedgelengths = [0.02,0.005];
    % Number of steps for curvature smoothing
    su2meshparam.curvatureSmoothing = 10;
end

function [strctout]=RSVS3D_GetOutputFiles(outDir, cellRecognize)
    % Gathers from a directory the outptu files specified by cell Recognize
    if ~exist('cellRecognize','var')
        cellRecognize = {'volu','^VoluMesh_.*.msh';
            'snaking','^SnakeMesh_.*.msh';
            'snake','^Snake_.*.3snk'};
    end
    filestruct = struct('name', '', 'path', '');
    strctout = struct();
    for ii = 1: size(cellRecognize,1)
        [pathStr, fileStr] = FindDirRegex(outDir,cellRecognize{ii,2},0);
        strctout.(cellRecognize{ii,1}) = filestruct;
        if numel(pathStr)==1
            strctout.(cellRecognize{ii,1}).path = pathStr{1};
            strctout.(cellRecognize{ii,1}).file = fileStr{1};
        elseif numel(pathStr)>1
            strctout.(cellRecognize{ii,1}).path = pathStr{1};
            strctout.(cellRecognize{ii,1}).file = fileStr{1};
            warning("Multiple (%i) files found matching : %s in:\n%s",...
                numel(pathStr) ,cellRecognize{ii,2}, outDir);
        else
            strctout.(cellRecognize{ii,1}).path = '';
            strctout.(cellRecognize{ii,1}).file = '';
        end
    end
end

function [meshFile] = SU2Mesh_RSVS3D(param3D, meshFolder, su2meshparam, ...
    meshName)
    
    if ~exist('su2meshparam','var')
        [su2meshparam] = RSVS3D_SU2Config();
    end
    if ~exist('meshName','var')
        meshName = 'mesh.su2';
    end
    
    jsonSu2Str = jsonencode(su2meshparam);
    su2ParamStr = [meshName,', ',jsonSu2Str];
    [param3D]=RSVS3D_CallInDirectory(param3D, meshFolder);
    [diroutfiles]=RSVS3D_GetOutputFiles(meshFolder);
    exportConf = {{'su2', su2ParamStr}};
    param3D = RSVS3D_CallExport(param3D, exportConf,diroutfiles.volu.path, ...
        diroutfiles.snaking.path, diroutfiles.snake.path);
    
    [outRSVS3D]=RSVS3D_Run(param3D);
    
    [meshPath, ~]= FindDirRegex(outRSVS3D.folder,...
        regexprep(meshName, '<pattern>','.*'),0);
    meshFile = meshPath{1};
end






















