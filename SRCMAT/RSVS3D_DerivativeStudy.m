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
        figs = repmat(struct('handle',[],'type','','data',[],'stats',[],'links',[]),[0,1]);
        
    end
    properties(Access=private)
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
    
    % function [resout]=RSVS3D_DerivativeStudy(entryPoint, varargin)
    %
    %     switch entryPoint
    %         case 'parse'
    %             [resout] = ParseFiles(varargin{:});
    %         case 'stats'
    %             [resout]=StatisticsDerivatives(varargin{:});
    %         case 'plot'
    %         otherwise
    %             error('Unrecognised entry point.')
    %     end
    % end
    methods
%% Interfaces
        function obj = RSVS3D_DerivativeStudy(derivDirectory, derivFileName)
            obj.derivDirectory = derivDirectory;
            obj.derivFileName = derivFileName;
        end
        
        function Parse(obj)
            [obj.resout] = obj.ParseFiles(obj.derivDirectory, ...
                obj.derivFileName);
            obj.MakeStale();
            if isempty(obj.layoutFields)
                
            end
        end
        
        function Stats(obj)
            [obj.statsout]=StatisticsDerivatives(obj.resout, obj.derivDirectory,...
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
%% Internals
    %% Parser
    % See ParseFiles.m
   %% Statisitics
   % see StatisticsDerivatives.m
   %% Plotting
        function []=PlotOptions(obj)
            
            obj.figs(end+1).type = plotType;
            switch plotType
                case 'sparsity'
                    % Spy the hessians and look at the matrix structures
                    PlotSparsity(obj, indices, figindex);
                case 'matsort'
                    % look at all the elements of a matrix sorted
                    % Distributions?
                case 'stat3d'
                    % 3D 'surface' for a statistic
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
        
        function []=PlotStat3D(obj, statList ,figindex)
            % Plots surfaces associated with the 
            if ~exist('figindex','var')
                figindex = numel(obj.figs);
                if figindex == 0
                    figindex = 1;
                end
            end
            
            obj.figs(figindex).stat = statList;
            
            nAx = numel(statList);
            obj.figs(figindex).handle = figure;
            
            for ii = 1:nAx
                ax(ii) = subplot(1,nAx,ii);
                hold on;
            end
            
            nExtraDim = numel(obj.layoutCoeffs)-2;
            layerThickness = cumprod(cellfun(@numel,obj.layoutCoeffs(3:end)));
            layerSize  = cellfun(@numel,obj.layoutCoeffs(3:end));
            nLayers = layerThickness(end);
            layerThickness = [1, layerThickness];
            [X,Y] = meshgrid(obj.layoutCoeffs{1:2});
            layout = obj.getLayout();
            colors = get(gca,'ColorOrder');
            colorcycle = @(x) colors(mod(x-1,size(colors,1))+1,:);
            
            escape_= @(str)regexprep(str,'_','\\_');
            for ii = 1:nAx
                tempArray = nan(size(obj.resout));
                currStat = statList{ii};
                if ~iscell(currStat)
                    currStat = regexp(currStat,'\.','split');
                end
                for jj = 1:numel(tempArray)
                    tempField = obj.statsout(jj);
                    
                    for kk = 1:numel(currStat)
                        tempField = tempField.(currStat{kk});
                    end
                    if numel(tempField)~=1 || ~isnumeric(tempField)
                        error(['The requested data did not reduce to a '...
                            'single number. Cannot compute surface.'])
                    end
                    tempArray(jj) = tempField;
                end
                varName = '';
                for kk = 1:numel(currStat)
                    varName= [varName , currStat{kk},'.'];
                end
                varName = escape_(varName(1:end-1));
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
            end
            TightenAxes(ax,'auto');
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
            obj.layoutFields=varargin;
            obj.staleLayout = true;
        end
        
        function reIndex(obj)
            cellInd = num2cell(1:numel(obj.resout));
            [obj.resout.index] = deal(cellInd{:});
            [obj.statsout.index] = deal(cellInd{:});
            obj.MakeStale();
        end
        
        function  layout = getLayout(obj)
            if obj.staleLayout
                if isempty(obj.layoutFields)
                    errstruct.message = ['layoutFields not set please set ' , ...
                        'using setLayoutFields(...)'];
                    errstruct.identifier = 'RSVS3D:DerivativeStudy:getLayout:unsetFields';
                    error(errstruct)
                end
                sizeLayout = zeros(1,3);
                flatLayout = cell(size(obj.layoutFields));
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
    end % methods
end % classdef














