function [statsout]=Structure2Statistics(resin, folderIn, derivFileName)
    % Function for the processing of standard metric 
    if ~exist('folderIn','var')
        folderIn = '';
        isSave = false;
    else
        isSave = true;
    end
    if ~exist('derivFileName','var')
        derivFileName = 'derivatives';
    end
    [statsout]=AnalyseAnyFormat(resin);
    if isSave
        save([folderIn,filesep,'stats_',derivFileName,'.mat'],'statsout')
    end
end% function


function [statsstruct]=AnalyseAnyFormat(derivIn)
    
    if isnumeric(derivIn)
        if numel(derivIn)>1
            statsstruct = ProcessNumeric(derivIn);
        else
            statsstruct = derivIn;
        end
    elseif isstruct(derivIn)
        %statsstruct = struct();
        if numel(derivIn)>0
            statsstruct = AnalyseScalarStructure(derivIn(1));
            statsstruct = repmat(statsstruct,[numel(derivIn),1]);
            for jj=2:numel(derivIn)
                statsstruct(jj) = AnalyseScalarStructure(derivIn(jj));
            end
        end
    elseif iscell(derivIn)
        statsstruct = cell(size(derivIn));
        for ii = 1:numel(derivIn)
            [statsstruct{ii}]=AnalyseAnyFormat(derivIn{ii});
        end
    end
    
    
end

function [statsstruct]=AnalyseScalarStructure(derivIn)
    statsstruct = struct();
    fields = fieldnames(derivIn)';
    for ii = 1:numel(fields)
        field = fields{ii};
        [statsstruct.(field)]=AnalyseAnyFormat(derivIn.(field));
    end
end % function

function [stats]=ProcessNumeric(arr, basefieldname)
    % Captures the numeric data in the array and reduces it to some basic
    % stats
    if ~exist('basefieldname','var')
        basefieldname= '';
    end
    siz = size(arr);
    dim = numel(siz);
    calcs = {'min','max','std', 'mean','median','std_no0', 'mean_no0',...
        'std_no0abs', 'mean_no0abs','norm'};
    if dim==2 && siz(1)==siz(2)
        calcs{end+1} = 'eig';
        calcs{end+1} = 'rcond';
    end
    
    stats = struct();
    for ii = 1:numel(calcs)
        calc = calcs{ii};
        switch calc
            case 'min'
                stats.([basefieldname,calc]) = min(arr(:));
            case 'max'
                stats.([basefieldname,calc]) = max(arr(:));
            case 'std'
                stats.([basefieldname,calc]) = std(arr(:));
            case 'mean'
                stats.([basefieldname,calc]) = mean(arr(:));
            case 'median'
                stats.([basefieldname,calc]) = median(arr(:));
            case 'std_no0'
                stats.([basefieldname,calc]) = std(arr(arr(:)~=0));
            case 'mean_no0'
                stats.([basefieldname,calc]) = mean(arr(arr(:)~=0));
            case 'std_no0abs'
                stats.([basefieldname,calc]) = std(abs(arr(arr(:)~=0)));
            case 'mean_no0abs'
                stats.([basefieldname,calc]) =  mean(abs(arr(arr(:)~=0)));
            case 'norm'
                stats.([basefieldname,calc]) = norm(arr);
            case 'rcond'
                stats.([basefieldname,calc]) = rcond(arr);
            case 'eig'
                eigs = eig(arr);
                stats.([basefieldname,calc]) = eigs;
                stats2 = ProcessNumeric(eigs, [calc,'_']);
                fields = fieldnames(stats2);
                for jj = 1:numel(fields)
                    field = fields{jj};
                    stats.(field) = stats2.(field);
                end
        end % switch
    end % for
end % function