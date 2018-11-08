

function [] = include_3DCheckGrid()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end

%% Base Grid operation

function [grid]=GridStructInit(varargin)
    
    grid=struct('vert',struct([]),'edge',struct([]),'surf',struct([]),...
        'volu',struct([]));
    grid.vert=struct('index',[],'coord',zeros([0 2]),'edgeind',zeros([1 0]));
    grid.edge=struct('index',[],'vertind',zeros([0 2]),'surfind',zeros([1 0]));
    grid.surf=struct('index',[],'fill',[],'edgeind',zeros([1 0]),'voluind',...
        zeros([0 2]));
    grid.volu=struct('index',[],'fill',[],'surfind',zeros([1 0]));
    
    if nargin==0
        grid.vert=repmat(grid.vert,[1,0]);
        grid.edge=repmat(grid.edge,[1,0]);
        grid.surf=repmat(grid.surf,[1,0]);
        grid.volu=repmat(grid.volu,[1,0]);
    elseif nargin==1
        grid.vert=repmat(grid.vert,[1,varargin{1}(1)]);
        grid.edge=repmat(grid.edge,[1,varargin{1}(2)]);
        grid.surf=repmat(grid.surf,[1,varargin{1}(3)]);
        grid.volu=repmat(grid.volu,[1,varargin{1}(4)]);
    elseif nargin==4
        grid.vert=repmat(grid.vert,[1,varargin{1}]);
        grid.edge=repmat(grid.edge,[1,varargin{2}]);
        grid.surf=repmat(grid.surf,[1,varargin{3}]);
        grid.volu=repmat(grid.volu,[1,varargin{4}]);
    else
        error('Unrecognised number of arguments.')
    end
    
end

function [grid]=LoadGridFromFile(filePath)
    [grid]=GridStructInit();
    fid=fopen(filePath,'r');
    
    str=fgetl(fid);
    nums=str2num(str);
    grid.vert(1).index=1;
    grid.vert=repmat(grid.vert,[1,nums(1)]);
    for ii=1:numel(grid.vert)
        str=fgetl(fid);
        nums=str2num(str);
        [grid.vert(ii)]=ParseVert(grid.vert(ii),nums);
    end
    
    
    str=fgetl(fid);
    nums=str2num(str);
    grid.edge(1).index=1;
    grid.edge=repmat(grid.edge,[1,nums(1)]);
    for ii=1:numel(grid.edge)
        str=fgetl(fid);
        nums=str2num(str);
        [grid.edge(ii)]=ParseEdge(grid.edge(ii),nums);
    end
    
    
    str=fgetl(fid);
    nums=str2num(str);
    grid.surf(1).index=1;
    grid.surf=repmat(grid.surf,[1,nums(1)]);
    for ii=1:numel(grid.surf)
        str=fgetl(fid);
        nums=str2num(str);
        [grid.surf(ii)]=ParseSurf(grid.surf(ii),nums);
    end
    
    
    str=fgetl(fid);
    nums=str2num(str);
    grid.volu(1).index=1;
    grid.volu=repmat(grid.volu,[1,nums(1)]);
    for ii=1:numel(grid.volu)
        str=fgetl(fid);
        nums=str2num(str);
        [grid.volu(ii)]=ParseVolu(grid.volu(ii),nums);
    end
    fclose(fid);
end

function [vert]=ParseVert(vert,nums)
    vert.index=nums(1);
    s=3;
    vert.edgeind=nums(s+(1:nums(s)));
    s=3+nums(s)+1;
    vert.coord=nums(s+(1:nums(s)));
end


function [vert]=ParseEdge(vert,nums)
    vert.index=nums(1);
    s=3;
    vert.vertind=nums(s+(1:nums(s)));
    s=3+nums(s)+1;
    vert.surfind=nums(s+(1:nums(s)));
end


function [vert]=ParseSurf(vert,nums)
    vert.index=nums(1);
    vert.fill=nums(2);
    s=6;
    vert.voluind=nums(s+(1:nums(s)));
    s=s+nums(s)+1;
    vert.edgeind=nums(s+(1:nums(s)));
end


function [vert]=ParseVolu(vert,nums)
    vert.index=nums(1);
    
    vert.fill=nums(2);
    s=6;
    vert.surfind=nums(s+(1:nums(s)));
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

function [voluOnRight]=FindRightSurface2D(grid)
    % Does not do anything at the moment simply returns one for each
    % surface. Tecplot survives without
    % To be updated if necessary.
    
    voluOnRight=ones(size(grid.edge));
    
    
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

function []=Check3Dgrid(grid,textLevel, plotNorm)
    
    if ~exist('textLevel','var');textLevel='ind';end
    if ~exist('plotNorm','var');plotNorm=0;end
    
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
    if plotNorm
        figure
        PlotSurfNormals(grid,axes)
    end 
    
    
    
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
            if ~isempty(cellStr{ii})
                text(mean(coord(arr,1)),mean(coord(arr,2)),mean(coord(arr,3)),cellStr{ii})
            end
        end
    else
        
        for ii=1:numel(grid.edge)
            arr=((ii-1)*2)+1:ii*2;
            plot(coord(arr,1),coord(arr,2),'o-')
            if ~isempty(cellStr{ii})
                text(mean(coord(arr,1)),mean(coord(arr,2)),cellStr{ii})
            end
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
            b = sort(b);
            if currVertSub(b(end))==currVertSub(end)
                b(end+1) = numel(currVertSub)-1;
            else
                b(end+1) = numel(currVertSub);
            end
            currVertSub=currVertSub(b);
            if currVertSub(end)~=currVertSub(1)
                currVertSub(1:2)=currVertSub(2:-1:1);
            end
            coord=vertcat(grid.vert(currVertSub).coord);
            mCoord=mean(coord,1);
            coord=(coord-repmat(mean(coord,1),[size(coord,1),1]))*0.95...
                +repmat(mCoord,[size(coord,1),1]);
            plot3(coord([1:end,1],1),coord([1:end,1],2),coord([1:end,1],3),...
                [vertexSymbols{mod(ii,numel(vertexSymbols))+1},'-']);
            if ~isempty(cellStr{ii})
                text(mCoord(:,1),mCoord(:,2),mCoord(:,3),cellStr{ii})
            end
        end
    else
        
        for ii=1:numel(grid.edge)
            arr=((ii-1)*2)+1:ii*2;
            plot(coord(arr,1),coord(arr,2),'o-')
            if ~isempty(cellStr{ii})
                text(mean(coord(arr,1)),mean(coord(arr,2)),cellStr{ii})
            end
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
            if ~isempty(cellStr{ii})
                text(mCoord(:,1),mCoord(:,2),mCoord(:,3),cellStr{ii})
            end
        end
    else
        
        for ii=1:numel(grid.edge)
            arr=((ii-1)*2)+1:ii*2;
            plot(coord(arr,1),coord(arr,2),'o-')
            if ~isempty(cellStr{ii})
                text(mean(coord(arr,1)),mean(coord(arr,2)),cellStr{ii})
            end
        end
    end
end

%% Surface normal

% Plot a normal to the face colored based on connectivity
% plot a centre connection between the face and the centre of the volume
% based on mean positions.

function []=PlotSurfNormals(grid, ax)
    vertInd=[grid.vert(:).index];
    edgeInd=[grid.edge(:).index];
    surfInd=[grid.surf(:).index];
    voluInd=[grid.volu(:).index];
    nVolu = numel(grid.volu);
    
    voluLines = cell(nVolu, 1);
    [voluLines{:}] = deal(zeros(0,3));
    voluColors=NiceColors();
    nColors = size(voluColors, 1);
    voluMCoord=zeros(nVolu, 3);
    for ii=1:nVolu
        voluMCoord(ii, :)=MeanVolume(grid.volu(ii), grid, edgeInd, vertInd,...
            surfInd);
    end
    
    for ii=1:numel(grid.surf)
        
        [normal, mCoord]=SurfNormal(grid.surf(ii), grid, edgeInd, vertInd);
        voluSurf = FindObjNum([],grid.surf(ii).voluind,voluInd);
        kk=-1;
        for jj = 1:numel(voluSurf)
            if voluSurf(jj)>0
                voluLines{voluSurf(jj)} = [ voluLines{voluSurf(jj)};
                    mCoord;mCoord+kk*normal;mCoord*nan;
                    mCoord;voluMCoord(voluSurf(jj),:);mCoord*nan];
            end
            kk=kk+2;
        end
        
    end
   hold on
    for ii=1:nVolu
        plot3(ax,voluLines{ii}(:,1),voluLines{ii}(:,2),voluLines{ii}(:,3),...
            '-', 'color', voluColors(mod(ii-1,nColors)+1,:));
    end
    
    
    
end

function [mCoord]=MeanVolume(volu, grid, edgeInd, vertInd, surfInd)
    currVertSub=FindObjNum([],...
        [grid.edge(...
        FindObjNum([],[grid.surf(...
        FindObjNum([],[volu.surfind],surfInd)).edgeind],edgeInd)...
        ).vertind],vertInd);
    [currVertSub,~,~]=unique(currVertSub);
    coord=vertcat(grid.vert(currVertSub).coord);
    mCoord=mean(coord,1);
end

function [normal, mCoord]=SurfNormal(surf, grid, edgeInd, vertInd)
    
    % Get vertices in right order
    currVertSub=FindObjNum([],...
                [grid.edge(...
                FindObjNum([],[surf.edgeind],edgeInd)...
                ).vertind],vertInd);
    [~,b,~]=unique(currVertSub);
    b = sort(b);
    if currVertSub(b(end))==currVertSub(end)
        b(end+1) = numel(currVertSub)-1;
    else
        b(end+1) = numel(currVertSub);
    end
    currVertSub=currVertSub(b);
    if currVertSub(end)~=currVertSub(1)
        currVertSub(1:2)=currVertSub(2:-1:1);
    end
    coords=vertcat(grid.vert(currVertSub).coord);
    coords(end+1,:)=coords(1,:);
    mCoord=mean(coords,1);
    
    normal = zeros(1,3);
    for ii=1:numel(currVertSub)
        normal = normal + TriPoint_Normal(mCoord,coords(ii,:),coords(ii+1,:));
    end
    
end

function [norm]=TriPoint_Normal(p0,p1,p2)
    
    norm = cross((p1-p0),(p2-p0));
    
end


function [rgb]=NiceColors()
    hexCol= {
    '#1f77b4'
    '#ff7f0e'
    '#2ca02c'
    '#d62728'
    '#9467bd'
    '#8c564b'
    '#e377c2'
    '#7f7f7f'
    '#bcbd22'
    '#17becf'};
    [ rgb ] = hex2rgb(hexCol);
end

function [ rgb ] = hex2rgb(hex,range)
% hex2rgb converts hex color values to rgb arrays on the range 0 to 1. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% SYNTAX:
% rgb = hex2rgb(hex) returns rgb color values in an n x 3 array. Values are
%                    scaled from 0 to 1 by default. 
%                    
% rgb = hex2rgb(hex,256) returns RGB values scaled from 0 to 255. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% EXAMPLES: 
% 
% myrgbvalue = hex2rgb('#334D66')
%    = 0.2000    0.3020    0.4000
% 
% 
% myrgbvalue = hex2rgb('334D66')  % <-the # sign is optional 
%    = 0.2000    0.3020    0.4000
% 
%
% myRGBvalue = hex2rgb('#334D66',256)
%    = 51    77   102
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myrgbvalues = hex2rgb(myhexvalues)
%    =   0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%        0.8000    0.6000    0.2000
%        0.2000    0.2000    0.9020
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myRGBvalues = hex2rgb(myhexvalues,256)
%    =   51    77   102
%       128   153   179
%       204   153    51
%        51    51   230
% 
% HexValsAsACharacterArray = {'#334D66';'#8099B3';'#CC9933';'#3333E6'}; 
% rgbvals = hex2rgb(HexValsAsACharacterArray)
% 
% * * * * * * * * * * * * * * * * * * * * 
% Chad A. Greene, April 2014
%
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% the improvement tips. In this update, the documentation now shows that
% the range may be set to 256. This is more intuitive than the previous
% style, which scaled values from 0 to 255 with range set to 255.  Now you
% can enter 256 or 255 for the range, and the answer will be the same--rgb
% values scaled from 0 to 255. Function now also accepts character arrays
% as input. 
% 
% * * * * * * * * * * * * * * * * * * * * 
% See also rgb2hex, dec2hex, hex2num, and ColorSpec. 
% 
%% Input checks:
assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.') 
if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
%% Tweak inputs if necessary: 
if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    
    % In case cell array elements are separated by a comma instead of a
    % semicolon, reshape hex:
    if isrow(hex)
        hex = hex'; 
    end
    
    % If input is cell, convert to matrix: 
    hex = cell2mat(hex);
end
if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end
if nargin == 1
    range = 1; 
end
%% Convert from hex to rgb: 
switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
    
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
end















