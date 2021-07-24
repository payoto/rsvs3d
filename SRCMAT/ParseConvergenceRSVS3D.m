function [resout] = ParseConvergenceRSVS3D(derivDirectory, derivFileName)
    [derivPaths, ~]=FindDir(derivDirectory, derivFileName, 0);
    resout = repmat(struct(),[0,1]);
    for ii = 1:numel(derivPaths)
        if ~strcmp(derivPaths{ii}(end-3:end),'.mat')
            [resout]=ReadDerivativeIn(derivPaths{ii}, resout);
        end
    end
    fileOut = [derivDirectory,filesep,'parsed_',derivFileName,'.mat'];
    save(fileOut,'resout');
end % function

function [resout]=ReadDerivativeIn(fileIn, resout)
    recordMarker = 'step';
    fid = fopen(fileIn, 'r');
    while ~feof(fid)
        str = fgetl(fid);
        if (numel(str)==0)
            continue
        end
        if ~isempty(regexp(str,recordMarker,'once'))
            % New structure
            resout(end+1).index = numel(resout)+1;
            varNames = regexp(str,',*','split','once');
            resout(end).time = str2num(varNames{2});
        else
            % new variable
            varNames = regexp(str,',*','split','once');

            % Read single number
            varTemp = deblank(regexprep(varNames{1},'^>\s*(.*):','$1'));
            var = matlab.lang.makeValidName(varTemp);
            resout(end).(var) = str2num(varNames{2}); %#ok<*ST2NM>

            if numel(resout(end).(var))==5 && varTemp(end)~='s'
                cellDat = {'norm','mean','std','max','min'};
                cellDat(2,:) = num2cell(reshape(resout(end).(var),[1,5]));
                resout(end).(var) = struct(cellDat{:});
            end

        end
    end
    fclose(fid);
end % function
