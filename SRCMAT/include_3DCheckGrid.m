

function [] = include_3DCheckGrid()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end

%% Grid Connectivity methods


function [vertInSurf]=FindVertInSurf(grid)
    vertInSurf=repmat(struct('vertind',[]),size(grid.surf));
    edgeInd=[grid.edge(:).index];
    vertInd=[grid.vert.index];
    
    for ii=1:numel(grid.surf)
        currVertSub=FindObjNum([],...
            [grid.edge(...
            FindObjNum([],[grid.surf(ii).edgeind],edgeInd)...
            ).vertind],vertInd);
        [~,b,~]=unique(currVertSub);
        currVertSub=currVertSub(sort(b));
        vertInSurf(ii).vertind=[grid.vert(currVertSub).index];
    end
    
    
end

function [voluOnRight]=FindRightVolu3D(grid)
    % Does not do anything at the moment simply returns one for each
    % surface. Tecplot survives without
    % To be updated if necessary.
    
    voluOnRight=ones(size(grid.surf));

    
end

%% Tecplot Output
function [cellStr]=PrepareCharForTecplot(cellStr)
    % Trims tecplot data to make sure it does not pass the maximum
    % character length of 30000
    numMax=30000;
    numChar=cellfun(@numel,cellStr);
    nNewPos=floor(numChar/numMax);
    sumNewPos=cumsum([0,nNewPos(1:end-1)]);
    ind=(1:numel(cellStr))+sumNewPos;
    %cellStr2=cell([1,numel(cellStr)+sum(nNewPos)]);
    cellStr(ind)=cellStr;
    for ii=find(nNewPos)
        jj=ii+sumNewPos(ii);
        for kk=jj:jj+nNewPos(ii)-1
            indCut=regexp(cellStr{kk}(numMax+1:end),' ','once')+numMax;
            cellStr{kk+1}=cellStr{kk}(indCut:end);
            cellStr{kk}(indCut:end)='';
        end
        
    end
    cellStr=cellStr(~cellfun(@isempty,cellStr));
end

%% Plot Generated grid

function []=Check3Dgrid(grid,textLevel)
    
    if ~exist('textLevel','var');textLevel='ind';end
    
    h=figure('Name','check grid');
    ax(1)=subplot(2,2,1);
    hold on
    Plot3DVert(grid,textLevel)
    ax(2)=subplot(2,2,2);
    hold on
    Plot3DEdge(grid,textLevel)
    ax(3)=subplot(2,2,3);
    hold on
    Plot3DSurf(grid,textLevel)
    ax(4)=subplot(2,2,4);
    hold on
    Plot3DVolume(grid,textLevel)
    Link = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget',...
        'XLim', 'YLim', 'ZLim'});
    setappdata(h, 'StoreTheLink', Link);
    
    
    
end

function []=Plot3DVert(grid,textLevel)
    
    coord=vertcat(grid.vert(:).coord);
    ind={grid.vert(:).index};
    edgeInd={grid.vert(:).edgeind};
    
    
    cellStr=cell(size(ind));
    switch textLevel
        case 'none'
            
        case 'all'
            for ii=1:numel(cellStr)
                cellStr{ii}=[int2str(ind{ii}),' (',int2str(edgeInd{ii}),')'];
            end
            
        case 'ind'
            for ii=1:numel(cellStr)
                cellStr{ii}=[int2str(ind{ii})];
            end
    end
    if size(coord,2)==3
        plot3(coord(:,1),coord(:,2),coord(:,3),'o')
        text(coord(:,1),coord(:,2),coord(:,3),cellStr)
    else
        
        plot3(coord(:,1),coord(:,2),zeros(size(coord(:,1))),'o')
        text(coord(:,1),coord(:,2),zeros(size(coord(:,1))),cellStr)
    end
end

function []=Plot3DEdge(grid,textLevel)
    
    coord=vertcat(grid.vert(FindObjNum([],[grid.edge.vertind],[grid.vert.index])).coord);
    ind={grid.edge(:).index};
    vertInd={grid.edge(:).vertind};
    surfInd={grid.edge(:).surfind};
    
    cellStr=cell(size(ind));
    switch textLevel
        case 'none'
            
        case 'all'
            for ii=1:numel(cellStr)
                cellStr{ii}=[int2str(ind{ii}),' (',int2str(vertInd{ii}),' ; ',...
                    int2str(surfInd{ii}),')'];
            end
            
        case 'ind'
            for ii=1:numel(cellStr)
                cellStr{ii}=[int2str(ind{ii})];
            end
            
    end
    if size(coord,2)==3
        for ii=1:numel(grid.edge)
            arr=((ii-1)*2)+1:ii*2;
            plot3(coord(arr,1),coord(arr,2),coord(arr,3),'o-')
            text(mean(coord(arr,1)),mean(coord(arr,2)),mean(coord(arr,3)),cellStr{ii})
        end
    else
        
        for ii=1:numel(grid.edge)
            arr=((ii-1)*2)+1:ii*2;
            plot(coord(arr,1),coord(arr,2),'o-')
            text(mean(coord(arr,1)),mean(coord(arr,2)),cellStr{ii})
        end
    end
end

function []=Plot3DSurf(grid,textLevel)
    
    edgeInd=[grid.edge(:).index];
    ind={grid.surf(:).index};
    vertInd=[grid.vert.index];
    edgeSurfInd={grid.surf(:).edgeind};
    voluSurfInd={grid.surf(:).voluind};
    
    cellStr=cell(size(ind));
    switch textLevel
        case 'none'
            
        case 'all'
            for ii=1:numel(cellStr)
                cellStr{ii}=[int2str(ind{ii}),' (',int2str(edgeSurfInd{ii}),' ; ',...
                    int2str(voluSurfInd{ii}),')'];
            end
            
        case 'ind'
            
            for ii=1:numel(cellStr)
                cellStr{ii}=[int2str(ind{ii})];
            end
    end
    vertexSymbols={'o','x','s','+','*','d','v'};
    if size(grid.vert(1).coord,2)==3
        for ii=1:numel(grid.surf)
            currVertSub=FindObjNum([],...
                [grid.edge(...
                FindObjNum([],[grid.surf(ii).edgeind],edgeInd)...
                ).vertind],vertInd);
            [~,b,~]=unique(currVertSub);
            currVertSub=currVertSub(sort(b));
            coord=vertcat(grid.vert(currVertSub).coord);
            mCoord=mean(coord,1);
            coord=(coord-repmat(mean(coord,1),[size(coord,1),1]))*0.95...
                +repmat(mCoord,[size(coord,1),1]);
            plot3(coord([1:end,1],1),coord([1:end,1],2),coord([1:end,1],3),...
                [vertexSymbols{mod(ii,numel(vertexSymbols))+1},'-']);
            text(mCoord(:,1),mCoord(:,2),mCoord(:,3),cellStr{ii})
        end
    else
        
        for ii=1:numel(grid.edge)
            arr=((ii-1)*2)+1:ii*2;
            plot(coord(arr,1),coord(arr,2),'o-')
            text(mean(coord(arr,1)),mean(coord(arr,2)),cellStr{ii})
        end
    end
end

function []=Plot3DVolume(grid,textLevel)
    edgeInd=[grid.edge(:).index];
    surfInd=[grid.surf(:).index];
    voluInd={grid.volu(:).index};
    surfVoluInd={grid.volu(:).surfind};
    vertInd=[grid.vert.index];
    edgeSurfInd={grid.surf(:).edgeind};
    voluFillInd={grid.volu(:).fill};
    
    cellStr=cell(size(voluInd));
    switch textLevel
        case 'none'
            
        case 'all'
            for ii=1:numel(cellStr)
                cellStr{ii}=[int2str(voluInd{ii}),' (',int2str(surfVoluInd{ii})...
                    ,' ; ',int2str(voluFillInd{ii}),')'];
            end
        case 'ind'
            for ii=1:numel(cellStr)
                cellStr{ii}=[int2str(voluInd{ii})];
            end
    end
    
    vertexSymbols={'o','x','s','+','*','d','v'};
    if size(grid.vert(1).coord,2)==3
        for ii=1:numel(grid.volu)
            currVertSub=FindObjNum([],...
                [grid.edge(...
                FindObjNum([],[grid.surf(...
                FindObjNum([],[grid.volu(ii).surfind],surfInd)).edgeind],edgeInd)...
                ).vertind],vertInd);
            [currVertSub,~,~]=unique(currVertSub);
            coord=vertcat(grid.vert(currVertSub).coord);
            mCoord=mean(coord,1);
            coord=(coord-repmat(mean(coord,1),[size(coord,1),1]))*0.95...
                +repmat(mCoord,[size(coord,1),1]);
            plot3([coord([1:end,1],1);mCoord(:,1)],[coord([1:end,1],2);...
                mCoord(:,2)],[coord([1:end,1],3);mCoord(:,3)],...
                [vertexSymbols{mod(ii,numel(vertexSymbols))+1}]);
            text(mCoord(:,1),mCoord(:,2),mCoord(:,3),cellStr{ii})
        end
    else
        
        for ii=1:numel(grid.edge)
            arr=((ii-1)*2)+1:ii*2;
            plot(coord(arr,1),coord(arr,2),'o-')
            text(mean(coord(arr,1)),mean(coord(arr,2)),cellStr{ii})
        end
    end
end

