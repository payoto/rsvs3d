function [varargout] = GetNestedStructureField(obj,varargin) %#ok<INUSL>
    
    varargout = cell(1,numel(varargin)*2);
    nVarArg=numel(varargin);
    nObj = numel(obj);
    
    strIndex = @(jj)['(',int2str(jj),')'];
    strRun = @(strEval,jj) ['varargout{ii}',strIndex(jj),' = obj',...
        strIndex(jj),strEval, ';'];
    for ii = 1:numel(varargin)
        % parse the characters
        if(ischar(varargin{ii}))
            varargin{ii} = regexp(varargin{ii},'[\.]','split');
        end
        strEval = '';
        nArg = numel(varargin{ii});
        if all(cellfun(@ischar,varargin{ii})) && all(isempty(regexp(varargin{ii},'[\(\)]','once')))

            for kk = 1:nArg
                strEval = [strEval,'.', varargin{ii}{kk}]; %#ok<*AGROW>
            end
            varargout{nVarArg+ii} = strEval;
        else % need to test for 
            for kk = 1:nArg
                if ischar(varargin{ii}{kk})
                    strEval = [strEval,'.', varargin{ii}{kk}];
                elseif isnumeric(varargin{ii}{kk}) && varargin{ii}{kk}>0 ...
                        && isfinite(varargin{ii}{kk})
                    strEval = [strEval,'(', int2str(varargin{ii}{kk}),')'];
                else
                    error('Unhandled type.')
                end
            end
        end
        varargout{nVarArg+ii} = ['obj',strEval];
        % Assign varargout
        jj=1;
        try
            eval(strRun(strEval, jj));
        catch ME
            ME2.message = [strRun(strEval, jj), '\n', ME.message];
            ME2.identifier = ME.identifier;
            ME2.stack = ME.stack;
            rethrow(ME2);
        end
        if nObj>1
            varargout{ii} = repmat(varargout{ii},[nObj,1]);
        end
        for jj = 2:nObj     
            try
                eval(strRun(strEval, jj));
            catch ME
                ME2.message = [strRun(strEval, jj), '\n', ME.message];
                ME2.identifier = ME.identifier;
                ME2.stack = ME.stack;
                rethrow(ME2);
            end
        end
        
    end
    
end