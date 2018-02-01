

function [] = include_CheckResultsLight()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end




function [movFrame]=CheckResults(iter,unstructured,oldGrid,snakposition,snaxel,makeMovie,volumefraction,borderVertices)
    global nDim domainBounds
    movFrame=[];
    if nDim==2
        figh=figure('Position',[100 100 1000 900]);
        axh=axes;
        hold on
        title(['Iteration  ',int2str(iter)],'fontsize',16)
        colString='bgcmyk';
        
        isEdgeSub=find(unstructured.edge.boundaryis1);
        for ii=1:length(isEdgeSub)
            PlotEdge(figh,axh,unstructured,isEdgeSub(ii),'bo')
        end
        
        isEdgeSub=find(unstructured.edge.boundaryis1);
        for ii=1:length(isEdgeSub)
            PlotEdge(figh,axh,unstructured,isEdgeSub(ii),'b-')
        end
        if exist('borderVertices','var')
            subVert=FindObjNum([],borderVertices,unstructured.vertex.index);
            for ii=1:length(subVert)
                PlotVert(figh,axh,unstructured,subVert(ii),'ro')
            end
        end
        isCellFull=find(unstructured.cell.fill);
        for ii=1:length( isCellFull)
            %PlotCell(figh,axh,unstructured, isCellFull(ii),'bs')
        end
        PlotSnaxel(figh,axh,snakposition,snaxel)
        %PlotSnaxelLoop(figh,axh,snakposition,snaxel)
        PlotSnaxelLoopDir(figh,axh,snakposition,snaxel)
        PlotSnaxelIndex(figh,axh,snakposition)
        
        if exist('volumefraction','var')
            oldCellIndUnstructInd=[oldGrid.cell(:).index];
            
            oldCellIndUnstructSub=FindObjNum(oldGrid.cell,oldCellIndUnstructInd);
            oldCellIndVolFracSub=FindObjNum(volumefraction,...
                oldCellIndUnstructInd,[volumefraction(:).oldCellInd]);
%             for ii=1:length(oldCellIndUnstructInd)
%                 
%                 coord=oldGrid.cell(oldCellIndUnstructSub(ii)).coord;
%                 % frac=volumefraction(oldCellIndVolFracSub(ii)).volumefraction...
%                 %     -volumefraction(oldCellIndVolFracSub(ii)).targetfill;
%                 frac=volumefraction(oldCellIndVolFracSub(ii)).oldCellInd;
%                 PlotVolFrac(figh,axh,coord,frac)
%             end
        end
        
        %         [normalcontourvec]=ContourNormal2(snaxel,snakposition);
        %         PlotContVec(figh,axh,snakposition,normalcontourvec)
        
        
        axis equal
        axis([domainBounds(1,1:2) domainBounds(2,1:2)])
        if makeMovie
            movFrame=getframe(figh);
        end
    end
    
end

function [figh]=CheckResultsLight(unstructured,snakposition,snaxel,figh)
    global nDim domainBounds
    
    if nDim==2
        if nargin==3
            figh=figure;
        else
            figure(figh)
        end
        axh=axes;
        hold on
        
        colString='bgcmyk';
        
        isEdgeIndex=find(unstructured.edge.boundaryis1);
        for ii=1:length(isEdgeIndex)
            %PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'bo')
        end
        
        isEdgeIndex=find(unstructured.edge.boundaryis0);
        for ii=1:length(isEdgeIndex)
            %PlotEdge(figh,axh,unstructured,isEdgeIndex(ii),'b-')
        end
        
        
        isCellFull=find(unstructured.cell.fill);
        for ii=1:length( isCellFull)
            %PlotCell(figh,axh,unstructured, isCellFull(ii),'bs')
        end
        %PlotSnaxel(figh,axh,snakposition)
        %PlotSnaxelLoop(figh,axh,snakposition,snaxel)
        PlotSnaxelLoopDir(figh,axh,snakposition,snaxel)
        %PlotSnaxelIndex(figh,axh,snakposition)
        
        %[normalcontourvec]=ContourNormal(snaxel,snakposition);
        %PlotContVec(figh,axh,snakposition,normalcontourvec)
        
        axis equal
        axis([domainBounds(1,1:2) domainBounds(2,1:2)])
    end
    
end

function []=PlotEdge(figh,axh,unstructured,subEdge,format)
    figure(figh)
    %axes(axh)
    
    vertices=unstructured.edge.vertexindex(subEdge,:);
    vertsub(1)=find(unstructured.vertex.index==vertices(1));
    vertsub(2)=find(unstructured.vertex.index==vertices(2));
    coord=unstructured.vertex.coord(vertsub,:);
    
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotVert(figh,axh,unstructured,subVert,format)
    figure(figh)
    %axes(axh)
    
    coord=unstructured.vertex.coord(subVert,:);
    
    plot(coord(:,1),coord(:,2),format)
    
end

function []=PlotSnaxel(figh,axh,snakposition,snaxel)
    % Plots the snaxels as arrows on the plot
    if numel(snaxel(1).v)==1
        for ii=1:length(snakposition)
            X(ii)=snakposition(ii).coord(1);
            Y(ii)=snakposition(ii).coord(2);
            U(ii)=snakposition(ii).vector(1)*snaxel(ii).v/40;
            V(ii)=snakposition(ii).vector(2)*snaxel(ii).v/40;
        end
    else
        for ii=1:length(snakposition)
            X(ii)=snakposition(ii).coord(1);
            Y(ii)=snakposition(ii).coord(2);
            U(ii)=snaxel(ii).v(1)/1000;
            V(ii)=snaxel(ii).v(2)/1000;
        end
    end
    figure(figh)
    axes(axh)
    quiver(X,Y,U,V,0,'r-')
    plot(X,Y,'rs')
    
end

function []=PlotContVec(figh,axh,snakposition,normalContVec)
    % Plots the snaxels as arrows on the plot
    snaxIndex=[snakposition(:).index];
    for ii=1:length(normalContVec)
        for jj=1:2
            workInd=normalContVec(ii).(['index',int2str(jj)]);
            %workInd=normalContVec(ii).vertex(jj);
            workSub(jj)=FindObjNum(snakposition,workInd,snaxIndex);
            coord(jj,1:2)=snakposition(workSub(jj)).coord;
        end
        coord=mean(coord);
        X(ii)=coord(1);
        Y(ii)=coord(2);
        U(ii)=normalContVec(ii).vector(1)/20;
        V(ii)=normalContVec(ii).vector(2)/20;
    end
    figure(figh)
    axes(axh)
    quiver(X,Y,U,V,0,'r-')
    
end

function []=PlotSnaxelIndex(figh,axh,snakposition)
    % Plots the snaxels as arrows on the plot
    figure(figh)
    axes(axh)
    for ii=1:length(snakposition)
        X(ii)=snakposition(ii).coord(1);
        Y(ii)=snakposition(ii).coord(2);
        U(ii)=snakposition(ii).vector(1)/5000;
        V(ii)=snakposition(ii).vector(2)/5000;
        str=num2str(snakposition(ii).index);
        text(X(ii)+U(ii),Y(ii)+V(ii),str)
    end
    
    
    
end

function []=PlotSnaxelLoop(figh,axh,snakposition,snaxel)
    % Plots the snaxels as arrows on the plot
    figure(figh)
    axes(axh)
    snaxInd=[snaxel(:).index];
    for jj=1:length(snaxel)
        line=[snaxel(jj).index,snaxel(jj).snaxnext];
        for ii=1:length(line)
            currSnaxSub=FindObjNum(snakposition,line(ii),snaxInd);
            X(ii)=snakposition(currSnaxSub).coord(1);
            Y(ii)=snakposition(currSnaxSub).coord(2);
        end
        plot(X,Y,'o--')
    end
    
end

function []=PlotSnaxelLoopDir(figh,axh,snakposition,snaxel)
    % Plots the snaxels as arrows on the plot
    figure(figh)
    axes(axh)
    snaxInd=[snaxel(:).index];
    for jj=1:length(snaxel)
        line=[snaxel(jj).index,snaxel(jj).snaxnext];
        for ii=1:length(line)
            currSnaxSub=FindObjNum(snakposition,line(ii),snaxInd);
            X(ii)=snakposition(currSnaxSub).coord(1);
            Y(ii)=snakposition(currSnaxSub).coord(2);
            
        end
        U=X(2)-X(1);
        
        V=Y(2)-Y(1);
        quiver(X(1),Y(1),U,V,0)
    end
    
end

function []=PlotVolFrac(figh,axh,coord,frac)
    figure(figh)
    axes(axh)
    if frac==0
        text(coord(:,1),coord(:,2),num2str(frac),'HorizontalAlignment','center')
    else
        text(coord(:,1),coord(:,2),num2str(frac,'%.1e'),'HorizontalAlignment','center')
    end
    hold on
end
