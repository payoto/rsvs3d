

function [] = include_SnakeParam()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end


%% From Snakes

function [unstructReshape]=ModifUnstructured(unstructured)
    % Reshapes the unstructureddata structure to b ein line with the shape
    % of "snakes"
    unstrucFields=fieldnames(unstructured);
    nFields=length(unstrucFields);
    
    for ii=1:nFields
        field1Fields=fieldnames(unstructured.(unstrucFields{ii}));
        nFields1=length(field1Fields);
        nObjects=length(unstructured.(unstrucFields{ii}).index);
        
        for jj=1:nObjects
            for kk=1:nFields1
                if ~isstruct(unstructured.(unstrucFields{ii}).(field1Fields{kk}))
                    
                    unstructReshape.(unstrucFields{ii})(jj).(field1Fields{kk})=...
                        unstructured.(unstrucFields{ii}).(field1Fields{kk})(jj,:);
                else
                    field2Fields=fieldnames(unstructured.(unstrucFields{ii}).(field1Fields{kk}));
                    nFields2=length(field2Fields);
                    
                    for ll=1:nFields2
                        unstructReshape.(unstrucFields{ii})(jj).(...
                            [field1Fields{kk},field2Fields{ll}])=...
                            unstructured.(unstrucFields{ii}).(...
                            field1Fields{kk}).(field2Fields{ll})(jj,:);
                    end
                end
            end
        end
    end
end

function [unstructured]=ModifReshape(unstructReshape)
    % Reshapes the unstructureddata structure to b ein line with the shape
    % of "snakes"
    unstrucFields=fieldnames(unstructReshape);
    nFields=length(unstrucFields);
    
    for ii=1:nFields
        field1Fields=fieldnames(unstructReshape.(unstrucFields{ii}));
        nFields1=length(field1Fields);
        nObjects=length(unstructReshape.(unstrucFields{ii}));
        
        
        for kk=1:nFields1
            unstructured.(unstrucFields{ii}).(field1Fields{kk})=...
                zeros([nObjects,length(unstructReshape.(unstrucFields{ii})(1).(field1Fields{kk}))]);
            for jj=1:nObjects
                %                 if ~isstruct(unstructReshape.(unstrucFields{ii}).(field1Fields{kk}))
                
                unstructured.(unstrucFields{ii}).(field1Fields{kk})(jj,:)...
                    =unstructReshape.(unstrucFields{ii})(jj).(field1Fields{kk});
                %                 else
                %                     field2Fields=fieldnames(unstructReshape.(unstrucFields{ii}).(field1Fields{kk}));
                %                     nFields2=length(field2Fields);
                %
                %                     for ll=1:nFields2
                %
                %                         unstructReshape.(unstrucFields{ii}).(...
                %                             field1Fields{kk}).(field2Fields{ll})(jj,:)= ...
                %                             unstructured.(unstrucFields{ii})(jj).(...
                %                             [field1Fields{kk},field2Fields{ll}]);
                %                     end
                %                 end
            end
        end
    end
end

function [leftMost]=LeftMostCorner(coord)
    % Returns the left most coordinate in a a set
    
    [xMin]=min(coord(:,1));
    iXMin=find(coord(:,1)==xMin);
    [~,iYMin]=min(coord(iXMin,2));
    leftMost=iXMin(iYMin);
    
end

function [isCCW]=CCWLoop(coord)
    % Checks if the order of points at the left most corner to determine the
    % direction of the loop.
    [mCoord,~]=size(coord);
    %coord(end-1:end,:)=[];
    isCCW=[];
    kk=0;
    while isempty(isCCW) && kk<4
        [leftMostCorner]=LeftMostCorner(coord);
        if ~isfinite(leftMostCorner)
            error('Left most corner is nan or inf')
        end
        switch leftMostCorner
            case 1
                precVert=mCoord;
                nextVert=leftMostCorner+1;
            case mCoord
                precVert=leftMostCorner-1;
                nextVert=1;
            otherwise
                precVert=leftMostCorner-1;
                nextVert=leftMostCorner+1;
        end

        precVec=coord(precVert,:)-coord(leftMostCorner,:);
        nextVec=coord(nextVert,:)-coord(leftMostCorner,:);
        precAngle=ExtractAngle360([-1 -1],precVec);
        nextAngle=ExtractAngle360([-1 -1],nextVec);


        if precAngle>nextAngle
            isCCW=true;
        elseif precAngle<nextAngle
            isCCW=false;
        else
            isCCW=[];
            kk=kk+1;
            coord=([0 1;-1 0]*coord')';
        end
    end
end

% Order BlockEdges might cause problems as the different versions were not
% consistant.
function [cellOrderedVertex,cellOrderedEdges]=...
        OrderBlockEdges(blockEdges,blockCellTrunc,edgeVertCellIntern)
    
    if ~exist('blockCellTrunc','var')
        blockCellTrunc=ones([length(blockEdges(:,1)),1]);
        edgeVertCellIntern=ones(1,4);
    end
    
    [mBE,~]=size(blockEdges);
    blockEdgesWorking=blockEdges;
    blockCellTruncWorking=blockCellTrunc;
    edgeList=1:mBE;
    
    % New array counters
    iCell=1;
    iEdge=0;
    % Old array locations
    ii=1;
    jj=1;
    while ~isempty(blockEdgesWorking)
        iEdge=iEdge+1;
        kk=abs(jj-3); % opposite column of jj
        % Save current Edge
        currentVertex=blockEdgesWorking(ii,jj);
        nextVertex=blockEdgesWorking(ii,kk);
        cellOrderedVertex{iCell}(iEdge,1)=currentVertex;
        cellOrderedVertex{iCell}(iEdge,2)=nextVertex;
        cellOrderedEdges{iCell}(iEdge)=edgeList(ii);
        
        % Delete current edge and edgeList entry from working set
        edgeList(ii)=[];
        blockEdgesWorking(ii,:)=[];
        currBlockCell=blockCellTruncWorking(ii);
        %Increment the counter variables
        
        blockCellTruncWorking(ii)=[];
        [ii,jj]=find(blockEdgesWorking==nextVertex);
        if length(ii)>1
            nAct=find(blockCellTruncWorking(ii)==currBlockCell);
            
            if isempty(nAct)
                tempBlock=blockCellTruncWorking(ii);
                nAct=find(~any(repmat(tempBlock,flip(size(tempBlock)))==...
                    repmat(tempBlock',size(tempBlock)) & ~eye(numel(ii))));
            end
            if numel(nAct)>1
                [currBlockCell]=ExploreEdgeCellConnectivity(nextVertex,...
                    currBlockCell,edgeVertCellIntern);
                nAct=find(blockCellTruncWorking(ii)==currBlockCell);
                if numel(nAct)>1
                   error('Failed to find a connected set of vertices')
                end
                
            end
            ii=ii(nAct);
            jj=jj(nAct);
            disp('Loops neighbouring at corner')
        end
        if isempty(ii) % reset loop if ii is not found
            % restart from the first unassigned edge
            ii=1;
            jj=1;
            % Increment the loop number
            iCell=iCell+1;
            % Restart teh edge count
            iEdge=0;
        elseif numel(ii)>1
            error('Initial Snake building process failed')
        end
    end
    
end

function [newTargCell]=ExploreEdgeCellConnectivity(nextVertex,targCell,...
        edgeVertCellIntern)
    % Should work with convex grids
    [ii,~]=find(edgeVertCellIntern(:,1:2)==nextVertex);
    edgeVertCellIntern=edgeVertCellIntern(ii,:);
    tempBlock=[targCell,reshape(edgeVertCellIntern(:,3:4)',[1,...
        numel(edgeVertCellIntern(:,3:4))])];
    connLog=repmat(tempBlock,flip(size(tempBlock)))==...
                    repmat(tempBlock',size(tempBlock)) & ~eye(numel(tempBlock));
    nCell=size(edgeVertCellIntern,1);
    connLogRed=zeros(nCell+1);
    connLogRed(1,2:end)=sum(reshape(connLog(1,2:end),[2,nCell]));
    connLogRed(2:end,1)=connLogRed(1,2:end)';
    for ii=1:size(edgeVertCellIntern,1)
        for jj=1:size(edgeVertCellIntern,1)
            connLogRed(ii+1,jj+1)=sum(sum(connLog(2+(ii-1)*2:2+(ii-1)*2+1,...
                2+(jj-1)*2:2+(jj-1)*2+1)));
        end
    end
    nIt=0;
    iiprec=1;
    ii=find(connLogRed(iiprec,:));
    while any(connLogRed(ii,[1:iiprec(end)-1,iiprec(end)+1:end])) && nIt<nCell
        nIt=nIt+1;
        iiprec(end+1)=ii;
        
        ii=find(connLogRed(ii,:));
        ii=ii(ii~=iiprec(end-1));
    end
    ii=ii-1;
    iiprec=iiprec(end)-1;
    jj=find(~any(connLog(max(2+(iiprec-1)*2:2+(iiprec-1)*2+1,1),max(2+(ii-1)*2:2+(ii-1)*2+1,1)),1));
    
    newTargCell=edgeVertCellIntern(ii,jj+2);
end

function [vecAngles]=ExtractAnglepm180(baseVector,testVector)
    % This function calculates the angle between vectors
    
    toComplex=[1;0+1i];
    baseAngle=angle(baseVector*toComplex);
    vecAngles=angle(testVector*toComplex)-baseAngle;
    vecAngles(vecAngles>pi)=vecAngles(vecAngles>pi)-2*pi;
    vecAngles(vecAngles<-pi)=vecAngles(vecAngles<-pi)+2*pi;
    
    
end

function [vecAngles]=ExtractAngle360(baseVector,testVector)
    % This function calculates the angle between vectors
    
    toComplex=[1;0+1i];
    baseAngle=angle(baseVector*toComplex);
    vecAngles=angle(testVector*toComplex)-baseAngle;
    vecAngles(vecAngles>(2*pi))=vecAngles(vecAngles>(2*pi))-2*pi;
    vecAngles(vecAngles<0)=vecAngles(vecAngles<0)+2*pi;
    
    
end

function [cellCentredGrid]=CellCentredGrid(refinedGrid)
    % Returns cell centred information about the grid being used
    
    cellCentredGrid=refinedGrid.cell;
    edgeCellInfo=vertcat(refinedGrid.edge(:).cellindex);
    
    origEdgeIndex=[refinedGrid.edge(:).index];
    origVertexIndex=[refinedGrid.vertex(:).index];
    origCellIndex=[refinedGrid.cell(:).index];
    
    for ii=1:length(cellCentredGrid)
        edgeCellLog=sum((edgeCellInfo==cellCentredGrid(ii).index),2)>0;
        totEdges=sum(edgeCellLog);
        cellCentredGrid(ii).edge(1:totEdges)=refinedGrid.edge(edgeCellLog);
        cellVertex=[cellCentredGrid(ii).edge(:).vertexindex];
        cellVertex=RemoveIdenticalEntries(cellVertex');
        cellVertexSub=FindObjNum(refinedGrid.vertex,cellVertex,origVertexIndex);
        totVertices=length(cellVertex);
        cellCentredGrid(ii).vertex(1:totVertices)=refinedGrid.vertex(cellVertexSub);
        
    end
end

function [vertexCentredGrid]=VertexCentredGrid(refinedGrid)
    % Returns cell centred information about the grid being used
    
    vertexCentredGrid=refinedGrid.vertex;
    vertexCellInfo=vertcat(refinedGrid.edge(:).vertexindex);
    
    origEdgeIndex=[refinedGrid.edge(:).index];
    origVertexIndex=[refinedGrid.vertex(:).index];
    origCellIndex=[refinedGrid.cell(:).index];
    
    for ii=1:length(vertexCentredGrid)
        vertCellLog=sum((vertexCellInfo==vertexCentredGrid(ii).index),2)>0;
        totEdges=sum(vertCellLog);
        vertexCentredGrid(ii).edge(1:totEdges)=refinedGrid.edge(vertCellLog);
        vertexCell=[vertexCentredGrid(ii).edge(:).cellindex];
        vertexCell=RemoveIdenticalEntries(vertexCell');
        vertexCell(vertexCell==0)=[];
        cellVertexSub=FindObjNum(refinedGrid.vertex,vertexCell,origCellIndex);
        totCells=length(vertexCell);
        vertexCentredGrid(ii).cell(1:totCells)=refinedGrid.cell(cellVertexSub);
        
    end
end

function [quotient,left]=IntegerQuotient(a,b)
    % Divides a by b and gives the integer result and the leftover
    % Works best for positive numbers
    
    quotient=floor(a/b);
    left=a-(floor(a/b)*b);
end

function surrogatePoints=PointGeneration(ranges,N_surpoints)
    % Produces an array of points containing N_surpoints in each dimension
    % combining every point with every dimesion
    %   RANGES: is a D*2 matrix containing the lower and upper bounds of each
    %           variable
    %        N: is the number of graduations in each dimension
    
    [m_ranges,~]=size(ranges);
    
    for ii=1:m_ranges
        if ranges(ii,1)~= ranges(ii,2)
            X_inter(:,ii)=linspace(ranges(ii,1),ranges(ii,2),N_surpoints);
        else
            X_inter(1:N_surpoints,ii)=ranges(ii,1);
        end
        
    end
    
    % Generation of points for RBF generation
    [Dim,~]=size(ranges);
    
    X_RBF=[];
    for ii=1:Dim
        [m_X,~]=size(X_RBF);
        m_X=max([m_X,1]);
        inter=[X_RBF,X_inter(1,ii)*ones(m_X,1)];
        for jj=2:N_surpoints
            
            % for each partial point already in corners this loop combines all
            % the values of the subsequent variable
            inter=[inter;X_RBF,X_inter(jj,ii)*ones(m_X,1)];
        end
        X_RBF=inter;
    end
    
    surrogatePoints=X_RBF;
    
end

%% Loop building

function [loop]=GenerateSnakStartLoop(gridrefined2,boundstr)
    
    isEdge=[gridrefined2.edge(:).(boundstr{1})];
    cond=boundstr{3};
    [loop]=OrderSurfaceVertexReshape(gridrefined2,isEdge,cond);
    
    [loop]=EdgeInCondForVertex(loop,gridrefined2,cond);
end

function [loopsnaxel]=ExtractSnaxelLoops(snaxel,param)
    
    [loopsnaxel]=OrderSurfaceSnaxel(snaxel);
    [loopsnaxel]=FinishLoops(loopsnaxel,param);
end

function [loopsnaxel]=FinishLoops(loopsnaxel,param)
    
    varExtract={'edgeFinish','TEShrink','LEShrink','axisRatio'};
    [edgeFinish,TEShrink,LEShrink,axisRatio]=ExtractVariables(varExtract,param);
    
    for ii=1:length(loopsnaxel)
        
        loopsnaxel(ii).snaxel.coord(:,2)=loopsnaxel(ii).snaxel.coord(:,2)*axisRatio;
    end
    
    switch edgeFinish
        case 'shrink'
            for ii=1:length(loopsnaxel)
                loopsnaxel(ii).snaxel.coord=ShrinkEdges(loopsnaxel(ii).snaxel.coord,LEShrink,TEShrink);
            end
        case 'sharpen'
            for ii=1:length(loopsnaxel)
                [loopsnaxel(ii).snaxel.coord]=SharpenEdges(loopsnaxel(ii).snaxel.coord,TEShrink,LEShrink);
            end
        case 'none'
            
    end
    
end

function [loopsnaxel]=OrderSurfaceSnaxel(snaxel,snaxPositions)
    % function extracting the snaxels into their separate loops
    if nargin==1
        global unstructglobal
        snaxPositions=PositionSnakes(snaxel,unstructglobal);
    end
    nSnax=length(snaxel);
    blockSegments=zeros(2*nSnax,2);
    for ii=1:nSnax
        for jj=0:1
            blockSegments(2*ii-jj,:)=[snaxel(ii).index,snaxel(ii).connectivity(jj+1)];
        end
    end
    
    cellSimilar=FindIdenticalVector(blockSegments);
    for ii=1:length(cellSimilar)
        blockEdgeIndex(ii)=cellSimilar{ii}(1);
    end
    blockEdges=blockSegments(blockEdgeIndex,:);
    % Order edges into closed loops
    [cellOrderedVertex]=OrderBlockEdges(blockEdges);
    snaxIndex=[snaxel(:).index];
    for ii=1:length(cellOrderedVertex)
        loopsnaxel(ii).snaxel.index=[cellOrderedVertex{ii}(:,1)];
        loopIndices=FindObjNum(snaxel,loopsnaxel(ii).snaxel.index,snaxIndex);
        loopsnaxel(ii).snaxel.coord=vertcat(snaxPositions(loopIndices).coord);
        %loopsnaxel(ii).edge.index=isEdgeIndex(cellOrderedEdges{ii});
    end
    
end

function [points]=ShrinkEdges(points,LEShrink,TEShrink)
    % Function which allows to make sharp trailing edges and leading edges
    [~,iLE]=min(points(:,1));
    [~,iTE]=max(points(:,1));
    
    [m,~]=size(points);
    
    iLEm1=iLE-1;
    [iLEm1]=IndexMod(iLEm1,m);
    iLEp1=iLE+1;
    [iLEp1]=IndexMod(iLEp1,m);
    
    iTEm1=iTE-1;
    [iTEm1]=IndexMod(iTEm1,m);
    iTEp1=iTE+1;
    [iTEp1]=IndexMod(iTEp1,m);
    
    points(iLEm1,2)=points(iLEm1,2)-LEShrink;
    points(iLEp1,2)=points(iLEp1,2)+LEShrink;
    
    points(iTEm1,2)=points(iTEm1,2)+TEShrink;
    points(iTEp1,2)=points(iTEp1,2)-TEShrink;
    
    if points(iLEm1,2)<=points(iLE+1,2)
        points([iLEm1,iLE+1],:)=[];
    end
    if points(iTEm1,2)>=points(iTEp1,2)
        points([iTEm1,iTEp1],:)=[];
    end
    
end

function [points]=SharpenEdges(points,TEShrink,LEShrink)
    % Function which allows to make sharp trailing edges and leading edges
    [~,iLE]=min(points(:,1));
    [~,iTE]=max(points(:,1));
    
    [m,~]=size(points);
    
    iLEm1=iLE-1;
    [iLEm1]=IndexMod(iLEm1,m);
    iLEp1=iLE+1;
    [iLEp1]=IndexMod(iLEp1,m);
    iLEm2=iLE-2;
    [iLEm2]=IndexMod(iLEm2,m);
    iLEp2=iLE+2;
    [iLEp2]=IndexMod(iLEp2,m);
    
    iTEm1=iTE-1;
    [iTEm1]=IndexMod(iTEm1,m);
    iTEp1=iTE+1;
    [iTEp1]=IndexMod(iTEp1,m);
    iTEm2=iTE-2;
    [iTEm2]=IndexMod(iTEm2,m);
    iTEp2=iTE+2;
    [iTEp2]=IndexMod(iTEp2,m);
    
    if LEShrink
        [points(iLEm1,:)]=Align3points(points([iLE,iLEm2],:),points(iLEm1,:));
        [points(iLEp1,:)]=Align3points(points([iLE,iLEp2],:),points(iLEp1,:));
    end
    if TEShrink
        [points(iTEm1,:)]=Align3points(points([iTE,iTEm2],:),points(iTEm1,:));
        [points(iTEp1,:)]=Align3points(points([iTE,iTEp2],:),points(iTEp1,:));
    end
    
    %     if points(iLEm1,2)<=points(iLE+1,2)
    %          points([iLEm1,iLE+1],:)=[];
    %     end
    %     if points(iTEm1,2)>=points(iTEp1,2)
    %          points([iTEm1,iTEp1],:)=[];
    %     end
    
end

function [pointAlign]=Align3points(line,point)
    
    pointAlign=point;
    
    pointAlign(2)=(line(1,2)-line(2,2))/(line(1,1)-line(2,1))...
        *(point(1)-line(2,1))+line(2,2);
    
    if isnan(pointAlign(2)) || ~isfinite(pointAlign(2))
        pointAlign=point;
    end
    
    
end

function [indMod]=IndexMod(ind,m)
    
    indMod=mod(ind-1,m)+1;
    
end

function [snakposition]=PositionSnakes(snaxel,unstructured)
    % Returns an array with Snaxel coordinates preceded by snaxel indices
    if numel(unstructured.vertex)==1
        vertIndex=unstructured.vertex.index;
        vertCoord=unstructured.vertex.coord;
    else
        vertIndex=[unstructured.vertex.index]';
        vertCoord=vertcat(unstructured.vertex.coord);
        
    end
    fromVertex=[snaxel(:).fromvertex];
    toVertex=[snaxel(:).tovertex];
    
    nSnaxel=length(snaxel);
    
    for ii=nSnaxel:-1:1
        iToVert=vertCoord(find(vertIndex==toVertex(ii)),:); %#ok<FNDSB> % extract vertex coordinates
        iFromVert=vertCoord(find(vertIndex==fromVertex(ii)),:); %#ok<FNDSB>
        
        snakposition(ii).index=snaxel(ii).index;
        snakposition(ii).coord=iFromVert+(iToVert-iFromVert)*snaxel(ii).d;
        snakposition(ii).vectornotnorm=(iToVert-iFromVert);
        snakposition(ii).vertInit=iFromVert;
        snakposition(ii).vector=(iToVert-iFromVert)/norm(iToVert-iFromVert);
    end
    
end

%% Grid Operation
% From ExecuteOptimisation

function [cellCentredCoarse]=CellCentredSnaxelInfo(snaxel,refinedGrid,...
        cellCentredFine,cellCentredCoarse,connecstruct)
    
%     oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
%                 {connecstruct.cell(:).newCellInd},...
%                 {connecstruct.cell(:).oldCellInd},'UniformOutput',false));
%     newIndsCell=[connecstruct.cell(:).newCellInd];   

    [oldIndsNewOrd,newIndsCell]=OldIndexToNewOrder(connecstruct.cell,'CellInd');
            
    
    newEdgeInd=[refinedGrid.edge(:).index];
    [cellCentredFine]=IdentifyCellSnaxel(snaxel,refinedGrid,cellCentredFine);
    
    % this line matches each Fine cell to its coarse cell in the
    % cellCentredGrid
    newToOldCell=FindObjNum([],oldIndsNewOrd(FindObjNum([],[cellCentredFine(:).index],...
        newIndsCell)),[cellCentredCoarse(:).index]);
    cellCentredCoarse(1).snaxel=struct([]);
    for ii=1:numel(newToOldCell)
        cellCentredCoarse(newToOldCell(ii)).snaxel=[cellCentredCoarse(newToOldCell(ii)).snaxel, ...
            cellCentredFine(ii).snaxel];
    end
    for ii=1:numel(cellCentredCoarse)
        cellCentredCoarse(ii).lSnax=0;
        cellCentredCoarse(ii).curvSnax=0;
        cellCentredCoarse(ii).curvNoEdge=0;
        cellCentredCoarse(ii).lSnaxNorm=0;
        coord=vertcat(cellCentredCoarse(ii).vertex(:).coord);
        
        cellCentredCoarse(ii).cellLength=max(coord)-min(coord);
        if ~isempty(cellCentredCoarse(ii).snaxel)
            [cellCentredCoarse(ii).snaxel]...
                =CompressSnaxelChain(cellCentredCoarse(ii).snaxel);
            % mark snaxel which are not on the border of the cell (ie who's edge has only one cell)
            for jj=1:numel(cellCentredCoarse(ii).snaxel)
                snaxCellNew=refinedGrid.edge((FindObjNum([],cellCentredCoarse(ii).snaxel(jj).edge,newEdgeInd))).cellindex;
                snaxCellOld=unique(oldIndsNewOrd(FindObjNum([],snaxCellNew,newIndsCell)));
                cellCentredCoarse(ii).snaxel(jj).isborder=numel(snaxCellOld)>1;
            end
        end
        
        

    end
    for ii=1:numel(cellCentredCoarse)
        if ~isempty(cellCentredCoarse(ii).snaxel)
            [cellCentredCoarse(ii).lSnax,cellCentredCoarse(ii).curvSnax,...
                cellCentredCoarse(ii).lSnaxNorm,cellCentredCoarse(ii).curvNoEdge]...
                =ExploreSnaxelChain(cellCentredCoarse(ii).snaxel,cellCentredCoarse(ii).cellLength);
        end
    end
    
end

function [lSnax,curvSnax,lSnaxNorm,curvSnaxNoBord]=ExploreSnaxelChain(snaxelpart,cellLength)
    
    [snaxInd]=[snaxelpart(:).index];
    [snaxPrec]=[snaxelpart(:).snaxprec];
    [snaxNext]=[snaxelpart(:).snaxnext];
    snaxCon=[snaxPrec;snaxNext];
    isVisited=false(size(snaxInd));
    lSnax=0;
    curvSnax=sum(sqrt(sum(vertcat(snaxelpart(:).curv).^2,2)));
    curvSnaxNoBord=sum(sqrt(sum(vertcat(snaxelpart(~logical([snaxelpart(:).isborder])).curv).^2,2)));
%     plotPoints= @(points) plot(points([1:end],1),points([1:end],2));
%     plotPoints(vertcat(snaxelpart(:).coord));
    while ~all(isVisited)
        
        startSub=find(~isVisited,1,'first');
        
        for isfwd=0:1
            currSub=startSub;
            flag=true;
            while flag
                isVisited(currSub)=true;
                nextSub=FindObjNum([],snaxCon(isfwd+1,currSub),snaxInd);
                if nextSub==0
                    break
                end
                coord1=snaxelpart(currSub).coord;
                coord2=snaxelpart(nextSub).coord;
                lSnax=lSnax+sqrt(sum((coord2-coord1).^2));
                lSnaxNorm=lSnax+sqrt(sum(((coord2-coord1)./cellLength).^2));
                currSub=nextSub;
                flag=~isVisited(currSub);
            end
        end
        
    end
    
end

function [cellCentredGrid]=IdentifyCellSnaxel(snaxel,refinedGrid,cellCentredGrid)
    % Extracts the snaxel data and matches it to the cells
    
    cellCentredGrid(1).snaxel=struct([]);
    [snakposition]=PositionSnakesStruct(snaxel,refinedGrid);
    [snaxel]=CalculateSnaxelCurvature(snaxel,snakposition);
    %     [snakPosInd]=ReferenceCompArray(snakposition,inf,'inf','index');
    %     [edgeInd]=ReferenceCompArray(refinedGrid.edge,inf,'inf','index');
    %     [cellInd]=ReferenceCompArray(cellCentredGrid,inf,'inf','index');
    %     [vertIndex]=ReferenceCompArray(refinedGrid.vertex,inf,'inf','index');
    %     [vertCoord]=ReferenceCompArrayVertCat(refinedGrid.vertex,inf,'inf','coord');
    
    snakPosInd=[snakposition(:).index];
    edgeInd=[refinedGrid.edge(:).index];
    cellInd=[cellCentredGrid(:).index];
    vertIndex=[refinedGrid.vertex(:).index];
    vertCoord=vertcat(refinedGrid.vertex(:).coord);
    
    
    
    %     snakPosFields='coord,vector,vectorprec,vectornext,normvector,edgelength';
    %     snakPosFields=VerticalStringArray(snakPosFields,',');
    snakPosFields={'coord','vector','vectorprec','vectornext','normvector','edgelength'};
    for ii=1:length(snaxel)
        snaxEdge=snaxel(ii).edge;
        snaxEdgeSub=FindObjNum(refinedGrid.edge,snaxEdge,edgeInd);
        snaxCells=refinedGrid.edge(snaxEdgeSub).cellindex;
        snaxCells=snaxCells(snaxCells~=0);
        snaxCellsSub=FindObjNum(cellCentredGrid,snaxCells,cellInd);
        snaxPosSub=FindObjNum(snakposition,snaxel(ii).index,snakPosInd);
        
        for jj=1:length(snaxCellsSub)
            iToVert=vertCoord(find(vertIndex==snaxel(ii).tovertex),:); %#ok<FNDSB> % extract vertex coordinates
            iFromVert=vertCoord(find(vertIndex==snaxel(ii).fromvertex),:); %#ok<FNDSB>
            
            snaxelCell=snaxel(ii);
            for kk=1:length(snakPosFields(:,1))
                fieldNam=deblank(snakPosFields{kk});
                snaxelCell.(fieldNam)=0;
            end
            
            snaxelCell.coord=snakposition(snaxPosSub).coord;
            snaxelCell.vector=snakposition(snaxPosSub).vector;
%             snaxelCell.vectorprec=snakposition(snaxPosSub).vectorprec;
%             snaxelCell.vectornext=snakposition(snaxPosSub).vectornext;
%             snaxelCell.normvector=snakposition(snaxPosSub).normvector;
            snaxelCell.edgelength=norm(iToVert-iFromVert);
            
            cellCentredGrid(snaxCellsSub(jj)).snaxel=...
                [cellCentredGrid(snaxCellsSub(jj)).snaxel,snaxelCell];
        end
    end
    
end

function [snakposition]=PositionSnakesStruct(snaxel,unstructured)
    % Returns an array with Snaxel coordinates preceded by snaxel indices
    vertIndex=[unstructured.vertex(:).index];
    vertCoord=vertcat(unstructured.vertex(:).coord);
    fromVertex=[snaxel(:).fromvertex];
    toVertex=[snaxel(:).tovertex];
    
    nSnaxel=length(snaxel);
    
    for ii=nSnaxel:-1:1
        iToVert=vertCoord(find(vertIndex==toVertex(ii)),:); %#ok<FNDSB> % extract vertex coordinates
        iFromVert=vertCoord(find(vertIndex==fromVertex(ii)),:); %#ok<FNDSB>
        
        snakposition(ii).index=snaxel(ii).index;
        snakposition(ii).coord=iFromVert+(iToVert-iFromVert)*snaxel(ii).d;
        snakposition(ii).vectornotnorm=(iToVert-iFromVert);
        snakposition(ii).vertInit=iFromVert;
        snakposition(ii).vector=(iToVert-iFromVert)/norm(iToVert-iFromVert);
    end
    
end

function [snaxel]=CalculateSnaxelCurvature(snaxel,snakposition)
    
    curvFunc=@(pi,pip1,pim1,s1,s2)(-pi*(s1+s2)+pip1*s2+pim1*s1)/(s1^2*s2+s2^2*s1);
    
    
    snaxInd=[snakposition(:).index];
    snaxPrecSub=FindObjNum([],[snaxel(:).snaxprec],snaxInd);
    snaxNextSub=FindObjNum([],[snaxel(:).snaxnext],snaxInd);
%     figure,hold on
    for ii=1:numel(snaxel)
        pi1=snakposition(ii).coord;
        pip1=snakposition(snaxNextSub(ii)).coord;
        pim1=snakposition(snaxPrecSub(ii)).coord;
        snaxel(ii).curv=curvFunc(pi1,pip1,pim1,norm(pi1-pip1),norm(pi1-pim1));
%         quiver(pi1(1),pi1(2),snaxel(ii).curv(1)/50,snaxel(ii).curv(2)/50)
%         plot(pi1(1),pi1(2),'r*')
    end
    
    
end

function [snaxel]=CalculateSnaxelTangent(snaxel,snakposition)
    
    if nargin==1
        snakposition=snaxel;
    end
    tanFunc=@(pi,pip1,pim1,s1,s2)(pi*(s1^2-s2^2)+pip1*s2^2-pim1*s1^2)/(s1^2*s2+s2^2*s1);
    
    
    snaxInd=[snakposition(:).index];
    snaxPrecSub=FindObjNum([],[snaxel(:).snaxprec],snaxInd);
    snaxNextSub=FindObjNum([],[snaxel(:).snaxnext],snaxInd);
    % figure,hold on
    for ii=1:numel(snaxel)
        pi1=snakposition(ii).coord;
        pip1=snakposition(snaxNextSub(ii)).coord;
        pim1=snakposition(snaxPrecSub(ii)).coord;
        snaxel(ii).tangent=tanFunc(pi1,pip1,pim1,norm(pi1-pip1),norm(pi1-pim1));
        snaxel(ii).tangent=snaxel(ii).tangent/sqrt(sum(snaxel(ii).tangent.^2));
%          quiver(pi1(1),pi1(2),snaxel(ii).tangent(1)/50,snaxel(ii).tangent(2)/50)
%          plot(pi1(1),pi1(2),'r*')
    end
    
    
end

function [snaxel]=CompressSnaxelChain(snaxel)
    
    snaxInd=[snaxel(:).index];
    
    [~,uniqSub]=unique(snaxInd);
    
    snaxel=snaxel(uniqSub);
    
end

function [oldIndsNewOrd,newInds]=OldIndexToNewOrder(connec,addstr,oldstrfull)
    % connec=oldGrid.connec.cell
    % ADDSTR is the suffix in the CONNEC structure to the old and new fields
    % [oldIndsNewOrd,newInds]=OldIndexToNewOrder(CONNEC,ADDSTR)
    % OLDSTRFULL allows to specify both the new and old fields in full
    % [oldIndsNewOrd,newInds]=OldIndexToNewOrder(CONNEC,NEWSTR,OLDSTR)
    
    newstr=['new'];
    oldstr=['old'];
    switch nargin
        case 2
            newstr=[newstr,addstr];
            oldstr=[oldstr,addstr];
        case 3
            newstr=addstr;
            oldstr=oldstrfull;
    end
        
    
    oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {connec(:).(newstr)},...
                {connec(:).(oldstr)},'UniformOutput',false));
            
    if nargout>1
        newInds=[connec(:).(newstr)];
    end
end

function [gridcoarsen,coarsenconnec]=CoarsenGrid(gridrefined,gridbase,gridconnec)
   % removes all the edges internal to an original grid 
   
   if nargin==1
       gridbase=gridrefined.base;
       gridconnec=gridrefined.connec;
       gridrefined=gridrefined.refined;
   end
   coarsenconnec=gridconnec;
   try
       [oldIndsNewOrd,newInds]=OldIndexToNewOrder(gridconnec.cell);
       [coarsenconnec.cell.new]=deal(coarsenconnec.cell.old);
       try
           
           [coarsenconnec.cell.newCellInd]=deal(coarsenconnec.cell.oldCellInd);
       catch
           
       end
   catch
       [oldIndsNewOrd,newInds]=OldIndexToNewOrder(gridconnec.cell,'CellInd');
       [coarsenconnec.cell.newCellInd]=deal(coarsenconnec.cell.oldCellInd);
    end
   
   gridcoarsen=gridrefined;
   gridcoarsen.cell=gridbase.cell;
   oldIndsNewOrd=[0,oldIndsNewOrd];
   sub=FindObjNum([],[gridcoarsen.edge.cellindex],[0,newInds]);
   rmEdge=false(size(gridcoarsen.edge));
   for ii=1:numel(gridcoarsen.edge)
       gridcoarsen.edge(ii).cellindex=oldIndsNewOrd(sub((2*(ii-1))+1:2*(ii)));
       rmEdge(ii)=gridcoarsen.edge(ii).cellindex(1)==gridcoarsen.edge(ii).cellindex(2);
   end
   
   gridcoarsen.edge=gridcoarsen.edge(~rmEdge);
   
end



