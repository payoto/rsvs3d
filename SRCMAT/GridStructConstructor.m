%#codegen

function [grid]=GridStructConstructor(dimGrid,outType)
    % This function generates a sample grid for 2D and 3D RSVS.
    % Naming convention: 4 letters + 3 letter prefix
    %close all
    
    grid=struct('vert',struct([]),'edge',struct([]),'surf',struct([]),...
        'volu',struct([]));
    
    grid.vert=struct('index',[],'coord',zeros([0 2]),'edgeind',zeros([1 0]));
    grid.edge=struct('index',[],'vertind',zeros([0 2]),'surfind',zeros([1 0]));
    grid.surf=struct('index',[],'fill',[],'edgeind',zeros([1 0]),'voluind',...
        zeros([0 2]));
    grid.volu=struct('index',[],'fill',[],'surfind',zeros([1 0]));
    
    switch numel(dimGrid)
        case 2
            error('Not coded yet')
            [grid]=BuildSquare(grid);
%     [cubegrid]=BuildCube(grid);
        case 3
            [grid]=BuildBlock(grid,dimGrid);
    end
    
    
    switch outType
        case 'fig'
            Check3Dgrid(grid)
        case 'tecplot'
            fileName=sprintf('TESTOUT\\TestGrid3D%s.plt',num2str(dimGrid,'%i_'));
            fid=fopen(fileName,'w');
            Grid2Tecplot(grid,fid);
        case ''
        case 'none'
        otherwise
            warning('Output Type not recognised')
    end
end

function [squaregrid]=BuildSquare(squaregrid)
    
    squaregrid.vert=repmat(squaregrid.vert,[1,4]);
    squaregrid.edge=repmat(squaregrid.edge,[1,4]);
    squaregrid.volu=repmat(squaregrid.volu,[1,0]);
    
    coord=[0 0; 1 0; 1 1; 0 1];
    connec=[1:4;mod(1:4,4)+1]';
    
    
    squaregrid.surf.index=1;
    squaregrid.surf.fill=0.5;
    squaregrid.surf.edgeind=[1:4];
    
    for ii=1:4
        
        squaregrid.vert(ii).coord=coord(ii,:);
        squaregrid.vert(ii).index=ii;
        squaregrid.vert(ii).edgeind=connec(mod(ii-2,4)+1,:);
        
        squaregrid.edge(ii).index=ii;
        squaregrid.edge(ii).vertind=connec(ii,:);
        squaregrid.edge(ii).surfind=[1 0];
        
        
    end
    
end

function [cubegrid]=BuildCube(cubegrid)
    
    [squaregrid]=BuildSquare(cubegrid);
    
    cubegrid.vert=repmat(cubegrid.vert,[1,8]);
    cubegrid.edge=repmat(cubegrid.edge,[1,12]);
    cubegrid.surf=repmat(cubegrid.surf,[1,6]);
    cubegrid.volu=repmat(cubegrid.volu,[1,1]);
    
    cubegrid.volu.index=1;
    cubegrid.volu.fill=0.5;
    cubegrid.volu.surfind=[1:6];
    % First Face z=0
    for ii=1:numel(squaregrid.vert)
        squaregrid.vert(ii).coord=[squaregrid.vert(ii).coord,0];
        squaregrid.edge(ii).surfind(2)=squaregrid.edge(ii).surfind(2)+ii+1;
        squaregrid.vert(ii).edgeind(end+1)=8+ii;
    end
    
    [cubegrid.vert(1:4).index]=deal(squaregrid.vert(:).index);
    [cubegrid.vert(1:4).coord]=deal(squaregrid.vert(:).coord);
    [cubegrid.vert(1:4).edgeind]=deal(squaregrid.vert(:).edgeind);
    
    [cubegrid.edge(1:4).index]=deal(squaregrid.edge(:).index);
    [cubegrid.edge(1:4).vertind]=deal(squaregrid.edge(:).vertind);
    [cubegrid.edge(1:4).surfind]=deal(squaregrid.edge(:).surfind);
    % Second face z=1
    for ii=1:numel(squaregrid.vert)
        squaregrid.vert(ii).coord(end)=1;
        squaregrid.edge(ii).surfind(1)=6;
        
        squaregrid.vert(ii).index=squaregrid.vert(ii).index+4;
        squaregrid.vert(ii).edgeind(1:2)=squaregrid.vert(ii).edgeind(1:2)+4;
        
        squaregrid.edge(ii).index=squaregrid.edge(ii).index+4;
        squaregrid.edge(ii).vertind=squaregrid.edge(ii).vertind+4;
    end
    
    [cubegrid.vert(5:8).index]=deal(squaregrid.vert(:).index);
    [cubegrid.vert(5:8).coord]=deal(squaregrid.vert(:).coord);
    [cubegrid.vert(5:8).edgeind]=deal(squaregrid.vert(:).edgeind);
    
    [cubegrid.edge(5:8).index]=deal(squaregrid.edge(:).index);
    [cubegrid.edge(5:8).vertind]=deal(squaregrid.edge(:).vertind);
    [cubegrid.edge(5:8).surfind]=deal(squaregrid.edge(:).surfind);
    % Connecting them up
    connec=[1:4;5:8]';
    for ii=1:4
        cubegrid.edge(ii+8).index=ii+8;
        cubegrid.edge(ii+8).vertind=connec(ii,:);
        cubegrid.edge(ii+8).surfind=[mod(ii-2,4)+2,ii+1];
    end
    
    % Surfaces
    squaregrid.surf.voluind=[1 0];
    squaregrid.surf.fill=[];
    cubegrid.surf(1)=squaregrid.surf;
    
    cubegrid.surf(6)=squaregrid.surf;
    cubegrid.surf(6).edgeind=cubegrid.surf(6).edgeind+4;
    cubegrid.surf(6).index=6;
    for ii=1:4
        cubegrid.surf(ii+1).index=ii+1;
        cubegrid.surf(ii+1).edgeind=[ii,mod(ii,4)+9,ii+4,ii+8];
        cubegrid.surf(ii+1).voluind=[1 0];
        
        
    end
    
end

function [cubegrid]=BuildBlock(cubegrid,dimGrid)
    % takes in dimGrid a vector specifying the number of cells in each
    % dimension [l m n] for (x y z) all growth is handled first along x
    % then y then z
    
    dimGrid=reshape(dimGrid,[1,numel(dimGrid)]);
    nVolu=prod(dimGrid);
    nVert=prod(dimGrid+1);
    surfProp=flip(eye(numel(dimGrid)));
    nSurf=reshape(prod(repmat(dimGrid,flip(size(dimGrid)))...
        +surfProp,2),size(dimGrid));
    edgeProp=(1-eye(numel(dimGrid)));
    nEdge=reshape(prod(repmat(dimGrid,flip(size(dimGrid)))...
        +edgeProp,2),size(dimGrid));
    
    cubegrid.vert=repmat(cubegrid.vert,[1,nVert]);
    cubegrid.edge=repmat(cubegrid.edge,[1,sum(nEdge)]);
    cubegrid.surf=repmat(cubegrid.surf,[1,sum(nSurf)]);
    cubegrid.volu=repmat(cubegrid.volu,[1,nVolu]);
    
    % Volumes
    [cubegrid]=BuildBlockVolumes(cubegrid,dimGrid,nSurf,surfProp);
    % Surfaces
    [cubegrid]=BuildBlockSurfaces(cubegrid,dimGrid,nSurf,surfProp,edgeProp,nEdge);
    % Edges
    [cubegrid]=BuildBlockEdges(cubegrid,dimGrid,nEdge,edgeProp,nSurf,surfProp);
    % Vertices
    [cubegrid]=BuildBlockVertices(cubegrid,dimGrid,edgeProp,nEdge);
end

function [cubegrid]=BuildBlockVolumes(cubegrid,dimGrid,nSurf,surfProp)
    
    
    incrSurf=cumsum([0,nSurf(1:2)]);
    
    matSurf=repmat(dimGrid,flip(size(dimGrid)))+surfProp;
    
    incPos=cumprod([ones([6,1]),[matSurf(:,1:2);matSurf(:,1:2)]],2);
    
    pos=zeros([1,3]);
    
    for ii=1:numel(cubegrid.volu)
        cubegrid.volu(ii).index=ii;
        cubegrid.volu(ii).fill=rand;
        % surfind dependent on position in each dimension
        pos(1,1)=mod(ii-1,dimGrid(1))+1;
        pos(1,2)=mod(floor((ii-1)/dimGrid(1)),dimGrid(2))+1;
        pos(1,3)=mod(floor((ii-1)/(dimGrid(1)*dimGrid(2))),dimGrid(3))+1;
        cubegrid.volu(ii).coord=(pos-0.5)./dimGrid;
        
        
        pos=pos-1;
        cubegrid.volu(ii).surfind=(sum([repmat(pos,[3,1]);repmat(pos,[3,1])+...
            surfProp].*incPos,2)+[incrSurf,incrSurf]')'+1;
        
    end
    
    %     Check Cell Ordering
%         figure, hold on
%         grid2.vert=cubegrid.volu;
%         [grid2.vert.coord]=deal(grid2.vert.coord);
%         [grid2.vert.edgeind]=deal(grid2.vert.surfind);
%         Plot3DVert(grid2);
    
    
end

function [cubegrid]=BuildBlockSurfaces(cubegrid,dimGrid,nSurf,surfProp,edgeProp,nEdge)
    % struct('index',[],'fill',[],'edgeind',zeros([1 0]),'voluind',...
    %    zeros([0 2]));
    
    
    matSurf=repmat(dimGrid,flip(size(dimGrid)))+surfProp;
    incrSurf=cumsum([0,nSurf(1:2)]);
    
    matEdge=repmat(dimGrid,flip(size(dimGrid)))+edgeProp;
    incrEdge=cumsum([0,nEdge(1:2)]);
    incPos=cumprod([ones([6,1]),[matEdge(:,1:2);matEdge(:,1:2)]],2);
    pos=zeros([1,3]);
    ind=1:3;
    for ii=1:numel(cubegrid.surf)
        cubegrid.surf(ii).index=ii;
        jplane=1+sum((ii)>cumsum(nSurf));
        dimGridCur=matSurf(jplane,:);
        pos(1)=mod(ii-1-incrSurf(jplane),dimGridCur(1))+1;
        pos(2)=mod(floor((ii-incrSurf(jplane)-1)/dimGridCur(1)),dimGridCur(2))+1;
        pos(3)=mod(floor((ii-incrSurf(jplane)-1)/(dimGridCur(1)*dimGridCur(2))),dimGridCur(3))+1;
        
        
         cubegrid.surf(ii).coord=(pos-(1+surfProp(jplane,:))*0.5)./dimGrid;
        % Volumeind assignement done
        boundaryFlag=~any(pos>dimGrid);
        cubegrid.surf(ii).voluind(1)=(((pos-1)*cumprod(...
            [1,dimGrid(1:2)])')+1)*boundaryFlag;
        pos2=pos-surfProp(jplane,:);
        boundaryFlag=~any(pos2<1);
        cubegrid.surf(ii).voluind(2)=(((pos2-1)*cumprod(...
            [1,dimGrid(1:2)])')+1)*boundaryFlag;
        
        pos=pos-1;
        mask=zeros([4,3]);
        currDim=ind(flip(ind)~=jplane);
        currDim2=[currDim,currDim];
        mask(2:3,currDim)=eye(2);
        cubegrid.surf(ii).edgeind=(sum((repmat(pos,[4,1])+mask).*...
            incPos(currDim2,:),2)+incrEdge(currDim2)')'+1;
        
%         surfLog=[repmat(pos,[3,1]);repmat(pos,[3,1])+mask];
%         surfLog=any((surfLog<0) | (surfLog>[matSurf;matSurf]-1) ,2)' | flip([1:3,1:3])==jplane;
%         cubegrid.edge(ii).surfind(surfLog)=[];
        
    end
    
    %     Check Surface Ordering
%         figure, hold on
%         grid2.vert=cubegrid.surf;
%         [grid2.vert.coord]=deal(grid2.vert.coord);
%         [grid2.vert.edgeind]=deal(grid2.vert.voluind);
%         Plot3DVert(grid2);
%         figure, hold on
%         Plot3DVert(grid2);
end

function [cubegrid]=BuildBlockEdges(cubegrid,dimGrid,nEdge,edgeProp,nSurf,surfProp)
    % struct('index',[],'vertind',zeros([0 2]),'surfind',zeros([1 0]));
    
    
    
    matEdge=repmat(dimGrid,flip(size(dimGrid)))+edgeProp;
    incrEdge=cumsum([0,nEdge(1:2)]);
    pos=zeros([1,3]);
    
    
    matSurf=repmat(dimGrid,flip(size(dimGrid)))+surfProp;
    incPos=cumprod([ones([6,1]),[matSurf(:,1:2);matSurf(:,1:2)]],2);
    incrSurf=cumsum([0,nSurf(1:2)]);
    
    ind=1:3;
    for ii=1:numel(cubegrid.edge)
        cubegrid.edge(ii).index=ii;
        jplane=sum(ii>incrEdge);
        dimGridCur=matEdge(jplane,:);
        pos(1)=mod(ii-1-incrEdge(jplane),dimGridCur(1))+1;
        pos(2)=mod(floor((ii-incrEdge(jplane)-1)/dimGridCur(1)),dimGridCur(2))+1;
        pos(3)=mod(floor((ii-incrEdge(jplane)-1)/(dimGridCur(1)*dimGridCur(2))),...
            dimGridCur(3))+1;
        
        
        % Volumeind assignement done
        cubegrid.edge(ii).vertind=((([pos;(pos+1-edgeProp(jplane,:))]-1)*cumprod(...
            [1,dimGrid(1:2)+1])')+1)';
        cubegrid.edge(ii).coord=(pos-1+(1-edgeProp(jplane,:))*0.5)./dimGrid;
        
        pos=pos-1;
        mask=zeros(3);
        mask(ind(flip(ind)~=jplane),ind(ind~=jplane))=-eye(2);
        cubegrid.edge(ii).surfind=(sum([repmat(pos,[3,1]);repmat(pos,[3,1])+...
            mask].*incPos,2)+[incrSurf,incrSurf]')'+1;
        
        surfLog=[repmat(pos,[3,1]);repmat(pos,[3,1])+mask];
        surfLog=any((surfLog<0) | (surfLog>[matSurf;matSurf]-1) ,2)' | flip([1:3,1:3])==jplane;
        cubegrid.edge(ii).surfind(surfLog)=[];
    end
    
    %     Check Surface Ordering
%          grid2.vert=cubegrid.edge;
%          [grid2.vert.edgeind]=deal(grid2.vert.surfind);
%          Plot3DVert(grid2);
%          figure, hold on
%          [grid2.vert.edgeind]=deal(grid2.vert.vertind);
%          Plot3DVert(grid2);
end

function [cubegrid]=BuildBlockVertices(cubegrid,dimGrid,edgeProp,nEdge)
    
    incrEdge=cumsum([0,nEdge(1:2)]);
    dimGridVert=dimGrid+1;
    matEdge=repmat(dimGrid,flip(size(dimGrid)))+edgeProp;
    edgeProp=1-edgeProp;
    
    incPos=cumprod([ones([6,1]),[matEdge(:,1:2);matEdge(:,1:2)]],2);
    %incPos=cumprod([1,dimGrid(1:2)]);
    
    for ii=1:numel(cubegrid.vert)
        cubegrid.vert(ii).index=ii;
        % surfind dependent on position in each dimension
        pos(1,1)=mod(ii-1,dimGridVert(1))+1;
        pos(1,2)=mod(floor((ii-1)/dimGridVert(1)),dimGridVert(2))+1;
        pos(1,3)=mod(floor((ii-1)/(dimGridVert(1)*dimGridVert(2))),dimGridVert(3))+1;
        %cubegrid.volu(ii).coord=pos-0.5;
        
        pos=pos-1;
        cubegrid.vert(ii).coord=pos./(dimGrid);
        cubegrid.vert(ii).edgeind=(sum([repmat(pos,[3,1]);repmat(pos,[3,1])-...
            edgeProp].*incPos,2)+[incrEdge,incrEdge]')'+1;
        
        edgeLog=[repmat(pos,[3,1]);repmat(pos,[3,1])-...
            edgeProp];
        edgeLog=any(edgeLog<0 | edgeLog>[matEdge;matEdge]-1,2)';
        cubegrid.vert(ii).edgeind(edgeLog)=[];
    end
    %figure, hold on
    %Plot3DVert(cubegrid);
end


