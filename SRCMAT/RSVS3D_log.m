classdef RSVS3D_log < RSVS3D_DerivativeStudy



    methods(Access=protected)

        function BuildStructure(obj,varargin)
            % Function which builds the data structure at the core of the
            % object.
            % This function should be overriden in subclasses to specialise
            % how the structure is built
            [obj.resout] = ParseConvergenceRSVS3D(obj.derivDirectory, ...
                obj.derivFileName);
        end
    end
    methods
        function obj = RSVS3D_log(derivDirectory, autoProcess)
            if ~exist('autoProcess','var')
                autoProcess = false;
            end
            obj.derivDirectory = derivDirectory;
            obj.derivFileName = 'convergence_.*.log';
            [~, fileNames]=FindDirRegex(derivDirectory, obj.derivFileName, 0);
            obj.derivFileName = regexprep(fileNames{1},'\..*','');
            if autoProcess
                obj.Parse();
                obj.Stats();
                obj.Plot('convergence');
            end
        end

        function [ax]=Plot(obj, plotType, varargin)
            ax=Plot@RSVS3D_DerivativeStudy(obj, plotType, varargin{:});
            switch plotType
                case 'convergence'
                    % Convergence
                    ax = obj.PlotConvergence();
                otherwise
                    return;
            end
            obj.figs(end).varargin = varargin;
        end

        function [ax]=PlotConvergence(obj)
            fig = figure('Name','convergence');
            ax = axes;
            hold on;
            colors = get(gca,'ColorOrder');
            colorcycle = @(x) colors(mod(x-1,size(colors,1))+1,:);
            markers = {'none','+','o','x','d'};
            markercycle = @(x) markers{mod(x-1,numel(markers))+1};
            steps = [obj.statsout.index];
            fields = {'constraintResidual', 'objectiveResidual' };
            for ii = 1:numel(fields)
                l(ii)=obj.EnvelopeLines(ax,fields{ii}, steps,...
                    colorcycle(ii), markercycle(ceil(ii/size(colors,1))));
            end
            legend(l);
            ax.YScale = 'log';
            ax.XLim = [min(steps), max(steps)];
            ax.YLim = [1e-20 1e10];
            print(fig,[obj.derivDirectory,filesep,obj.derivFileName,'.png'],...
                '-dpng','-r300');
        end

        function mainLine=EnvelopeLines(obj,ax,fieldName, x,color, marker)

            statList = {{fieldName,'.',{'mean','min','max'}}};
            statList  = RSVS3D_DerivativeStudy.ParseStatLists(statList);

            [means, mins, maxs]=GetNestedStructureField(obj.statsout, statList{:});

            l= plot(ax,x, means);
            mainLine = l;
            l(2)= plot(ax,x, mins);
            l(3)= plot(ax,x, maxs);
            mainLine.DisplayName = fieldName;
            mainLine.Marker = marker;
            [l.Color] = deal(color);
            [l(2:3).LineStyle] = deal('--');
        end
    end
end
