classdef RSVS3D_DerivativeStudy < handle
    
    properties
        % Raw data
        resout;
        % statistics data set
        statsout;
        % Directory
        derivDirectory;
        % derivFileName
        derivFileName;
        % Figure plots
        figs = repmat(struct('handle',[],'type','','data',[],'stats',[],'links',[],'varargin',{}),[0,1]);
        
    end
    properties %(Access=protected)
        % list of indices
        indexList = [];
        % is the list of indices stale?
        staleIndex = true;
        % Fields used for layouts
        layoutFields = {};
        % 3 D layout of coefficients
        layoutData = zeros(0,0,0);
        % lists of unique coefficients
        layoutCoeffs = cell(0,3);
        staleLayout = true;
    end
    
    methods(Access=protected)
        
        function BuildStructure(obj,varargin)
            % Function which builds the data structure at the core of the
            % object.
            % This function should be overriden in subclasses to specialise
            % how the structure is built
            [obj.resout] = ParseFiles(obj.derivDirectory, ...
                obj.derivFileName);
            for ii = 1 : numel(obj.resout)
                obj.resout(ii).HLag = obj.resout(ii).HObj ...
                    + obj.resout(ii).HConstr;
            end
        end
    end
    
    methods
%% Interfaces

        function obj = RSVS3D_DerivativeStudy(derivDirectory, derivFileName)
            if nargin==2
                obj.derivDirectory = derivDirectory;
                obj.derivFileName = derivFileName;
            else
                obj.derivDirectory = '';
                obj.derivFileName = '';
            end
        end
        
        
        function Parse(obj,varargin)
            % 
            obj.BuildStructure(varargin{:})
            obj.MakeStale();
            if isempty(obj.layoutFields)
                obj.setLayoutFields();
            end
        end
        
        function Stats(obj)
            [obj.statsout]=Structure2Statistics(obj.resout, obj.derivDirectory,...
                obj.derivFileName);
        end % function
        
        function Save(obj)
            fileOut = [obj.derivDirectory,filesep,'object_',obj.derivFileName,'.mat'];
            save(fileOut,'obj');
        end
        
        function SaveLight(obj)
            objlight = RSVS3D_DerivativeStudy(obj.derivDirectory, obj.derivFileName);
            objlight.figs = obj.figs;
            
            fileOut = [objlight.derivDirectory,filesep,'objectlight_',...
                objlight.derivFileName,'.mat'];
            save(fileOut,'objlight');
        end

   %% Plotting
        function [ax]=Plot(obj, plotType, varargin)
            ax = struct();
            obj.figs(end+1).type = plotType;
            switch plotType
                case 'sparsity'
                    % Spy the hessians and look at the matrix structures
                    obj.PlotSparsity(varargin{1}, numel(obj.figs));
                case 'matsort'
                    % look at all the elements of a matrix sorted
                    % or Distributions?
                case 'stat3d'
                    % 3D 'surface' for a statistic
                    ax=obj.PlotStat3D(varargin{:});
            end
            obj.figs(end).varargin = varargin; 
        end
        
        function RePlot(obj)
            for ii = 1:numel(obj.figs)
                obj.Plot(obj.figs(ii).type,obj.figs(ii).varargin{:});
            end
        end
        
        function ClearPlots(obj)
            if numel(obj.figs)>1
                obj.figs = repmat(obj.figs(1),[0,1]);
            end
        end
        
        function []=PlotSparsity(obj, indices, figindex)
            % plots the sparsity of a matrix
            if ~exist('figindex','var')
                figindex = numel(obj.figs);
                if figindex == 0
                    figindex = 1;
                end
            end
            
            [pos]=obj.recordpos(indices);
            obj.figs(figindex).data = indices;
            obj.figs(figindex).handle = figure('Name', ['Sparsity ',...
                int2str(indices)], 'position', [200,200,1100,600]);
            ax(1) = subplot(1,2,1);
            hold on;
            ax(2) = subplot(1,2,2);
            hold on;
            k = 1;
            for p = pos(:)'
                hlag = obj.resout(p).HConstr+obj.resout(p).HObj;
                s = size(hlag);
                r = symrcm(hlag);
                
                [i,j] = find(hlag);
                i = s(1)-i+1;
                l(k,1) = plot(ax(1), i , j,...
                    '.','markersize',5,...
                    'DisplayName',['H ', int2str(obj.resout(p).index)]);
                [i,j] = find(hlag(r,r));
                i = s(1)-i+1;
                l(k,2) = plot(ax(2), i , j,...
                    '.','markersize',5,...
                    'DisplayName',['H ', int2str(obj.resout(p).index)]);
                k=k+1;
            end
            legend(l(:,1))
            legend(l(:,2))
        end
        
        function [ax]=PlotStat3D(obj, statList, varargin)
            % Plots surfaces associated with the statistics
            p = inputParser;
            addRequired(p,'statList',@iscell);
            addOptional(p,'figindex',max(numel(obj.figs),1));
            addParameter(p,'axesSettings',{}, @(x) iscell(x) && ~mod(numel(x),2));
            addParameter(p,'surfSettings',{}, @(x) iscell(x) && ~mod(numel(x),2));
            parse(p,statList,varargin{:});
            inputVars = fieldnames(p.Results);
            for ii = 1:numel(inputVars)
                eval([inputVars{ii},'=p.Results.(inputVars{ii});']);
            end
            
            obj.figs(figindex).stat = statList;
            
            % PreProcessing of statList, these can be "packed" to have a recursive
            % pattern definition.
            statList = RSVS3D_DerivativeStudy.ParseStatLists(statList);
            
            nAx = numel(statList);
            obj.figs(figindex).handle = figure;
            % Build the plots
            for ii = 1:nAx
                ax(ii) = subplot(1,nAx,ii);
                hold on;
            end
            % Prepare calculations
            layout = obj.getLayout();
            nExtraDim = numel(obj.layoutCoeffs)-2;
            layerThickness = cumprod(cellfun(@numel,obj.layoutCoeffs(3:end)));
            layerSize  = cellfun(@numel,obj.layoutCoeffs(3:end));
            layerThickness = [1, layerThickness];
            nLayers = layerThickness(end);
            [X,Y] = meshgrid(obj.layoutCoeffs{1:2});
            colors = get(gca,'ColorOrder');
            colorcycle = @(x) colors(mod(x-1,size(colors,1))+1,:);
            
            escape_= @(str)regexprep(str,'_','\\_');
            for ii = 1:nAx
                % Concatenates a nested structure field
                [tempArray,varName] = GetNestedStructureField(...
                    obj.statsout,statList{ii});
                varName = escape_(varName);
                
                s=[];
                for jj = 1:nLayers
                    extraDimSub = num2cell(mod(floor(...
                        ones(1,nExtraDim)*(jj-1)./layerThickness(1:end-1)),...
                        layerSize)+1);
                    
                    Z = tempArray(layout(:,:,extraDimSub{:}));
                    dataSetName = '';
                    for kk = 1:numel(extraDimSub)
                        dataSetName = [dataSetName , obj.layoutFields{kk+2},...
                            ': ', num2str(obj.layoutCoeffs{kk+2}(extraDimSub{kk})),...
                            '; '];
                    end
                    dataSetName = dataSetName(1:end-2);
                    s(jj) = surf(ax(ii),X,Y,Z, 'DisplayName', dataSetName,...
                        'FaceColor', colorcycle(jj));
                end
                ax(ii).Title.String = ['field: ', varName];
                ax(ii).ZAxis.Label.String = varName;
                ax(ii).XAxis.Label.String = escape_(obj.layoutFields{1});
                ax(ii).YAxis.Label.String = escape_(obj.layoutFields{2});
                legend(s)
                for jj= 1:numel(s)
                    SetNestedStructureField(s(jj),surfSettings{:}); %#ok<IDISVAR,USENS>
                end
            end
            TightenAxes(ax,'auto');
            for ii= 1:numel(ax)
                SetNestedStructureField(ax(ii),axesSettings{:}); %#ok<IDISVAR,USENS>
            end
            
            obj.figs(figindex).links=linkprop(ax,{'XLim','YLim','CameraPosition'});
        end
        
        
    %% Getters Setters and accessors
        function [pos]=recordpos(obj,indices, recalc)
            if ~exist('recalc','var')
                recalc = false;
            end
            if obj.staleIndex || recalc
                obj.indexList = [obj.resout.index];
            end
            pos = FindObjNum([],indices, obj.indexList);
        end
        
        function MakeStale(obj)
            obj.staleIndex = true;
            obj.staleLayout= true; 
        end
        
        function setLayoutFields(obj,varargin)
            
            if numel(varargin)==0
                fields = fieldnames(obj.resout);
                for ii=1:numel(fields)
                    if numel(obj.resout(1).(fields{ii}))==1 ...
                            && ~strcmp('index',fields{ii})
                        varargin{end+1} = fields{ii};
                    end
                end
            end
            
            obj.layoutFields=varargin;
            obj.staleLayout = true;
        end
        
        function reIndex(obj)
            cellInd = num2cell(1:numel(obj.resout));
            [obj.resout.index] = deal(cellInd{:});
            [obj.statsout.index] = deal(cellInd{:});
            obj.MakeStale();
        end
        
        function layout = getLayout(obj)
            if obj.staleLayout
                if isempty(obj.layoutFields)
                    errstruct.message = ['layoutFields not set please set ' , ...
                        'using setLayoutFields(...)'];
                    errstruct.identifier = 'RSVS3D:DerivativeStudy:getLayout:unsetFields';
                    error(errstruct)
                end
                sizeLayout = zeros(1,3);
                flatLayout = cell(size(obj.layoutFields));
                obj.layoutCoeffs = cell(size(obj.layoutFields));
                for ii = 1:numel(obj.layoutFields)
                    [obj.layoutCoeffs{ii},~, flatLayout{ii}] = ...
                        unique([obj.resout(:).(obj.layoutFields{ii})]);
                    sizeLayout(ii) = numel(obj.layoutCoeffs{ii});
                end
                obj.layoutData = zeros(sizeLayout);
                obj.layoutData(sub2ind(sizeLayout,flatLayout{:})) = 1:numel(obj.resout);
                obj.staleLayout = false;
            end
            layout = obj.layoutData;
        end
        
        function obj_small = downSelect(obj, selectFunc, newName)
            % Selects only part of the structure to use
            if~exist('newName','var')
                 newName = 'part';
            end
            newDirBase = [obj.derivDirectory,filesep,newName];
            flag = true; 
            ii = 0;
            while flag
                ii = ii + 1;
                newDir = [newDirBase,'_',int2str(ii)];
                flag = exist(newDir, 'dir');
            end
            obj_small = RSVS3D_DerivativeStudy(newDir,...
                [obj.derivFileName,'_',newName,'_',int2str(ii)]);
            obj_small.resout = obj.resout(arrayfun(selectFunc, obj.resout));
            % only makes the directory after succesfull splitting.
            mkdir(newDir)
            
        end
    end
    methods(Static)
        function [statList] = ParseStatLists(statList)
            ii = 1;
            while ii<=numel(statList)
                if(iscell(statList{ii}))
                    nCell = numel(statList{ii});
                    if nCell>1
                        if ischar(statList{ii}{1}) && ischar(statList{ii}{2})
                            statList{ii}{1} = [statList{ii}{1},statList{ii}{2}];
                            statList{ii}(2) = [];
                        elseif ischar(statList{ii}{1})
                            temp = cell(numel(statList{ii}{2}),1);
                            for jj = 1:numel(statList{ii}{2})
                                temp{jj}=[statList{ii}{1},statList{ii}{2}{jj}];
                            end
                            statList{ii}{1}=temp;
                            statList{ii}(2) = [];
                        elseif ischar(statList{ii}{2})
                            temp = cell(numel(statList{ii}{1}),1);
                            for jj = 1:numel(statList{ii}{1})
                                temp{jj}=[statList{ii}{1}{jj},statList{ii}{2}];
                            end
                            statList{ii}{1}=temp;
                            statList{ii}(2) = [];
                        else
                            n1 = numel(statList{ii}{1});
                            n2 = numel(statList{ii}{2});
                            temp = cell(n1*n2,1);
                            for j1 = 1:n1
                            for j2 = 1:n2
                                jj = j1+(j2-1)*n1;
                                temp{jj}=[statList{ii}{1}{j1},statList{ii}{2}{j2}];
                            end
                            end
                            if(numel(temp)~=n1*n2)
                                error('Should not happen')
                            end
                            statList{ii}{1}=temp;
                            statList{ii}(2) = [];
                        end
                    else
                        nCell = numel(statList{ii}{1});
                        statList(end+1:end+nCell-1)= statList{ii}{1}(2:end);
                        statList{ii} = statList{ii}{1}{1};
                    end
                else
                    ii = ii+1;
                end
            end
        end

    end % methods
end % classdef














