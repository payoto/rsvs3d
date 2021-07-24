function [] = include_GridCheck()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.

    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)

end


function [figh]=CheckGrid(gridreshape,axRatio,figh,format,pText)
    %domainBounds=[0 1 ; 0 1];

    [unstructured]=ModifReshape(gridreshape);

    if nargin<2
        axRatio=1;
    end
    if nargin<3
        figh=figure;
        axh=axes;
    else
        axh=findobj(figh,'type','axes');
    end
    if nargin<4
        format='b-';
    end
    if nargin<5
        pText=true;
    end

    hold on


    isEdgeSub=find(unstructured.edge.index);
    for ii=1:length(isEdgeSub)
        PlotEdgeGrid(figh,axRatio,unstructured,isEdgeSub(ii),format,pText)
    end

%     isEdgeSub=find(~unstructured.edge.boundaryis0);
%     for ii=1:length(isEdgeSub)
%         PlotEdgeGrid(figh,axh,unstructured,isEdgeSub(ii),'b-')
%     end


    axis equal
    %axis([domainBounds(1,1:2) domainBounds(2,1:2)])


end

function []=PlotEdgeGrid(figh,axRatio,unstructured,subEdge,format,pText)
    figure(figh)
    %axes(axh)

    vertices=unstructured.edge.vertexindex(subEdge,:);
    vertsub(1)=find(unstructured.vertex.index==vertices(1));
    vertsub(2)=find(unstructured.vertex.index==vertices(2));
    coord=unstructured.vertex.coord(vertsub,:);
    coord(:,2)=coord(:,2)*axRatio;
    plot(coord(:,1),coord(:,2),format)
    if pText
    text(mean(coord(:,1)),mean(coord(:,2)),int2str(unstructured.edge.index(subEdge)),'color','b')
    text(mean(coord(1,1)),mean(coord(1,2)),int2str(vertices(1)),'color','g')
    text(mean(coord(2,1)),mean(coord(2,2)),int2str(vertices(2)),'color','g')
    end
end

function []=PlotLoopGrid(figh,axh,loop,indexLoop,format)
    figure(figh)
    axes(axh)


    coord=loop(indexLoop).vertex.coord;

    plot(coord(:,1),coord(:,2),format)

end

function []=PlotSubDivGrid(figh,axh,loop,indexLoop,format)
    figure(figh)
    axes(axh)


    coord=loop(indexLoop).subdivision;

    plot(coord(:,1),coord(:,2),format)

end

function []=PlotCellGrid(figh,axh,unstructured,indexCell,format)
    figure(figh)
    axes(axh)


    coord=unstructured.cell.coord(indexCell,:);

    plot(coord(:,1),coord(:,2),format)

end

function []=PlotCellCentredGrid(cellCentredGrid,format,figh)

    plotPoints= @(points,format) plot(points([1:end,1],1),points([1:end,1],2),format);
    if nargin<2
        format='+';
    end
    if nargin<3
        figh=figure;
    end
    figure(figh);
    hold on
    for ii=1:length(cellCentredGrid)


        coord=vertcat(cellCentredGrid(ii).vertex(:).coord);

        plotPoints(coord,format)
        text(mean(coord(:,1)),mean(coord(:,2)),int2str(cellCentredGrid(ii).index))
    end



end
