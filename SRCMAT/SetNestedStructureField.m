function [obj] = SetNestedStructureField(obj,varargin)
    
    if mod(numel(varargin),2)
        error('Field name-value pairs are expected. An even number of varargs should be provided');
    end
    
    for ii = 1:2:numel(varargin)-1
        if(ischar(varargin{ii}))
            varargin{ii} = regexp(varargin{ii},'[\.]','split');
        end
        if all(cellfun(@ischar,varargin{ii})) && all(isempty(regexp(varargin{ii},'[\(\)]','once')))
            nArg = numel(varargin{ii});
            if nArg==1
                obj.(varargin{ii}{1}) = varargin{ii+1};
            elseif nArg==2
                obj.(varargin{ii}{1}).(varargin{ii}{2}) = varargin{ii+1};
            elseif nArg==3
                obj.(varargin{ii}{1}).(varargin{ii}{2}).(varargin{ii}{3})...
                    = varargin{ii+1};
            else
                strEval = 'obj';
                for kk = 1:nArg
                    strEval = [strEval,'.', varargin{ii}{kk}];
                end
                strEval = [strEval, '= varargin{ii+1};'];
                try
                    eval(strEval);
                catch ME
                    ME.message = [strEval, '\n', ME.message];
                    rethrow(ME);
                end
            end
        else
            strEval = 'obj';
            nArg = numel(varargin{ii});
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
            strEval = [strEval, '= varargin{ii+1};'];
            try
                eval(strEval);
            catch ME
                ME.message = [strEval, '\n', ME.message];
                rethrow(ME);
            end
                
        end
    end
    
end