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
                    warning(['Cause of error : ', strEval])
                    rethrow(ME);
                end
            end
        else
            strEval = 'obj';

            if isnumeric(varargin{ii+1})
                openArray = '(';
                closeArray = ')';
            else
                openArray = '{';
                closeArray = '}';
            end

            nArg = numel(varargin{ii});
            prevIsNum = false;

            for kk = 1:nArg
                if ischar(varargin{ii}{kk})
                    if prevIsNum
                        strEval = [strEval,closeArray];
                    end
                    strEval = [strEval,'.', varargin{ii}{kk}];
                elseif isnumeric(varargin{ii}{kk}) && varargin{ii}{kk}>0 ...
                        && isfinite(varargin{ii}{kk})
                    if ~prevIsNum
                        strEval = [strEval,openArray];
                    else
                        strEval = [strEval,','];
                    end
                    strEval = [strEval, int2str(varargin{ii}{kk})];
                    prevIsNum = true;
                else
                    error('Unhandled type.')
                end
            end
            if prevIsNum
                strEval = [strEval,closeArray];
            end
            strEval = [strEval, '= varargin{ii+1};'];
            try
                eval(strEval);
            catch ME
                warning(['Cause of error : ', strEval])
                rethrow(ME);
            end

        end
    end

end
