function [statsout]=Structure2Statistics(resin, folderIn, derivFileName)
    %     resout = repmat(struct(),[0,1]);
    if ~exist('folderIn','var')
        folderIn = '';
        isSave = false;
    else
        isSave = true;
    end
    if ~exist('derivFileName','var')
        derivFileName = 'derivatives';
    end
    for ii=1:numel(resin)
        statsout(ii) = AnalyseDerivativeRecord(resin(ii));
    end
    if isSave
        save([folderIn,filesep,'stats_',derivFileName,'.mat'],'statsout')
    end
end% function

function [statsstruct]=AnalyseDerivativeRecord(derivIn)
    statsstruct = struct();
    fields = fieldnames(derivIn)';
    for ii = 1:numel(fields)
        field = fields{ii};
        if numel(derivIn.(field))>1
            statsstruct.(field) = ProcessNumeric(derivIn.(field));
        else
            statsstruct.(field) = derivIn.(field);
        end
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
    calcs = {'min','max','std', 'mean','std_no0', 'mean_no0',...
        'std_no0abs', 'mean_no0abs','norm'};
    if dim==2 && siz(1)==siz(2)
        calcs{end+1} = 'eig';
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