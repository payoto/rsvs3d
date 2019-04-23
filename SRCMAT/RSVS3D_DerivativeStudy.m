function [resout]=RSVS3D_DerivativeStudy(entryPoint, varargin)
    
    switch entryPoint
        case 'parse'
            [resout] = ParseFiles(derivDirectory, derivFileName);
        case 'stats'
        case 'plot'
        otherwise
            error('Unrecognised entry point.')
    end
end

function [resout] = ParseFiles(derivDirectory, derivFileName)
    [derivPaths, ~]=FindDir(derivDirectory, derivFileName, 0);
    resout = repmat(struct(),[0,1]);
    for ii = 1:numel(derivPaths)
        [resout]=ReadDerivativeIn(derivPaths{ii}, resout);
    end
    fileOut = [derivDirectory,filesep,'parsed_',derivFileName,'.mat'];
    save(fileOut,'resout');
end

function [resout]=ReadDerivativeIn(fileIn, resout)
    
    recordMarker = 'GenerateDerivatives';
    fid = fopen(fileIn, 'r');
    while ~feof(fid)
        str = fgetl(fid);
        if (numel(str)==0)
            continue
        end
        if ~isempty(regexp(str,recordMarker,'once'))
            % New structure
            resout(end+1).index = numel(resout)+1;
        elseif ~isempty(regexp(str(1),'[a-zA-Z]', 'once'))
            % new variable
            varNames = regexp(str,'[ ,]*','split');
            if isempty(regexp(str,',', 'once'))
                % Read single number
                resout(end).(varNames{1}) = str2double(varNames{2});
            else
                % Read array
                lengthArray = str2double(varNames{2});
                resout(end).(varNames{1}) = zeros(lengthArray,1);
                nRead = 0;
                nLines = 0;
                while nRead<lengthArray
                    str = fgetl(fid);
                    numDouble = numel(regexp(str,' '));
                    resout(end).(varNames{1})(nRead+1:nRead+numDouble) = ...
                        str2num(str); %#ok<ST2NM>
                    nRead = nRead + numDouble;
                    nLines = nLines + 1;
                end
                
                if (nRead==lengthArray)
                    resout(end).(varNames{1})= ...
                        reshape(resout(end).(varNames{1}),[numDouble,nLines])';
                else
                    warning('Could not reshape member %s', varNames{1})
                end
                
            end
        end
        
    end
    fclose(fid);
end