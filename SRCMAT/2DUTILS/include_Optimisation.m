function [] = include_Optimisation()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end


%%
function [constrVal]=NacaOuterLimit0012(gridrefined,paramoptim)
    
    [constrVal]=NacaOuterLimit4d(gridrefined,paramoptim,'0012');
    %{
    naca4t=@(x,t,c,xMin)  5*t*c*(0.2969*sqrt((x-xMin)/c)-0.1260*((x-xMin)/c)...
        -0.3516*((x-xMin)/c).^2+0.2843*((x-xMin)/c).^3-0.1036*((x-xMin)/c).^4);
    integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*...
        (tDistrib(1:end-1)+tDistrib(2:end))/2]);
    
    warning('Will only work with square or rectangular grids')
    varExtract={'axisRatio'};
    axisRatio=ExtractVariables(varExtract,paramoptim.parametrisation);
    
    
    t=0.12;
    cellCentredGrid=CellCentreGridInformation(gridrefined);
    isActive=logical([cellCentredGrid(:).isactive]);
    actVerts=[cellCentredGrid(isActive).vertex];
    coord=vertcat(actVerts(:).coord);
    
    coord(:,2)=coord(:,2)*axisRatio;%;
    
    xPos=RemoveIdenticalEntries(coord(:,1));
    xMax=max(xPos);
    xMin=min(xPos);
    fillSub=zeros([1,sum(isActive)]);
    reqFrac=zeros([1,sum(isActive)]);
    actCellSub=find(isActive);
    for ii=1:numel(actCellSub)
        
        cellCoords=vertcat(cellCentredGrid(actCellSub(ii)).vertex(:).coord);
        cellCoords(:,2)=cellCoords(:,2)*axisRatio;
        [cornerCoord]=IdentifyCorners(cellCoords);
        posMin=min(cornerCoord);
        posMax=max(cornerCoord);
        
        x=linspace(posMin(1),posMax(1),200);
        tDistrib=naca4t(x,t,(xMax-xMin),xMin);
        y=min(max(tDistrib,posMin(2)),posMax(2))-min(max(-tDistrib,posMin(2)),posMax(2));
        %plot(x,y+posMin(2))
        vol=integr(x,y);
        fillSub(ii)=ii;
        reqFrac(ii)=vol(end)/cellCentredGrid(actCellSub(ii)).volume/axisRatio;
    end
    constrVal={fillSub,reqFrac};
    %}
end

function [constrVal]=NacaOuterLimit4d(gridrefined,paramoptim,nacaStr)
    a4_open=0.1015;
    a4_closed=0.1036;
    naca4t=@(x,t,c,xMin,a4,teps)  5*t*c*(0.2969*sqrt((x-xMin)/c)-0.1260*((x-xMin)/c)...
        -0.3516*((x-xMin)/c).^2+0.2843*((x-xMin)/c).^3-a4*((x-xMin)/c).^4)+((x-xMin)/c)*teps;
    
    naca4c=@(x,m,p,c,xMin) [m/p^2*(2*p*((x((x-xMin)<(p*c))-xMin)/c)-((x((x-xMin)<(p*c))-xMin)/c).^2),...
        m/(1-p)^2*((1-2*p)+2*p*((x((x-xMin)>=(p*c))-xMin)/c)-((x((x-xMin)>=(p*c))-xMin)/c).^2)];
    %     ((x(x>=(p*c))-xMin)/c)
    %
    %     [m*x(x<(cp*c))/p^2.*(2*p-x(x<(p*c))) ;...
    %         m*(1-x(x>=(p*c))/c)/(1-p)^2.*(1+x(x>=(p*c))-2*p)];
    
    teps=5.48e-04/2/0.8;
    integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*...
        (tDistrib(1:end-1)+tDistrib(2:end))/2]);
    
    warning('Will only work with square or rectangular grids')
    varExtract={'axisRatio'};
    axisRatio=ExtractVariables(varExtract,paramoptim.parametrisation);
    [m,p,t,refFlag]=ReadNacaString(nacaStr);
    if m==0 || p==0
        p=0.0001;
        m=0;
    end
    
    cellCentredGrid=CellCentreGridInformation(gridrefined);
    isActive=logical([cellCentredGrid(:).isactive]);
    actVerts=[cellCentredGrid(isActive).vertex];
    coord=vertcat(actVerts(:).coord);
    %figure, hold on
    coord(:,2)=coord(:,2)*axisRatio;%;
    
    xPos=RemoveIdenticalEntries(coord(:,1));
    xMax=1; %max(xPos);
    xMin=0; %min(xPos);
    fillSub=zeros([1,sum(isActive)]);
    reqFrac=zeros([1,sum(isActive)]);
    actCellSub=find(isActive);
    for ii=1:numel(actCellSub)
        
        cellCoords=vertcat(cellCentredGrid(actCellSub(ii)).vertex(:).coord);
        cellCoords(:,2)=cellCoords(:,2)*axisRatio;
        [cornerCoord]=IdentifyCorners(cellCoords);
        posMin=min(cornerCoord);
        posMax=max(cornerCoord);
        
        %x=linspace(posMin(1),posMax(1),200);
        x=min(max(linspace(posMin(1),posMax(1),200),xMin),xMax);
        if min(x)==xMin
            m2X=min(x(x~=xMin))-xMin;
            if ~isempty(m2X)
                xPosMin=(x==xMin);
                x(xPosMin)=min((cumsum(xPosMin(xPosMin))-1)*m2X*3/sum(xPosMin)+xMin,xMax);
            end
        end
        if max(x)==xMax
            m2X=max(x(x~=xMax))-xMax;
            if ~isempty(m2X)
                xPosMin=(x==xMax);
                x(xPosMin)=max((cumsum(xPosMin(xPosMin))-1)*m2X/sum(xPosMin)+xMax,xMin);
            end
        end
        tDistrib=naca4t(x,t,(xMax-xMin),xMin,a4_closed,teps);
        cDistrib=naca4c(x,m,p,(xMax-xMin),xMin);
        y=min(max(cDistrib+tDistrib,posMin(2)),posMax(2))-min(max(cDistrib-tDistrib,posMin(2)),posMax(2));
        %plot(x,cDistrib+y+posMin(2))
        vol=integr(x,y);
        fillSub(ii)=ii;
        reqFrac(ii)=vol(end)/cellCentredGrid(actCellSub(ii)).volume/axisRatio;
    end
    constrVal={fillSub,reqFrac};
end

function [nacaLoops]=NacaMultiTopo(nPtsPloop,nacaStr)
    
    % xTop and xBot need to be normalised
    a4_open=0.1015;
    a4_closed=0.1036;
    
    naca45t=@(x,t,c,xMin,a4,teps)  5*t*c*(0.2969*sqrt((x-xMin)/c)-0.1260*((x-xMin)/c)...
        -0.3516*((x-xMin)/c).^2+0.2843*((x-xMin)/c).^3-a4*((x-xMin)/c).^4)+((x-xMin)/c)*teps;
    naca4c=@(x,m,p,c,xMin) [m/p^2*(2*p*((x((x-xMin)<(p*c))-xMin)/c)-((x((x-xMin)<(p*c))-xMin)/c).^2),...
        m/(1-p)^2*((1-2*p)+2*p*((x((x-xMin)>=(p*c))-xMin)/c)-((x((x-xMin)>=(p*c))-xMin)/c).^2)];
    teps=5.48e-04/2/0.8; % true @ corner=1e-5
    
    % Nine digit and more if for multi element airfoil
    % uses separator ";" between foils
    % uses _ separator between parameters
    % first digit is the number of airfoils
    % then each airfoil is placed as follows:
    % [NACA number l (in 1/10ths of chord) alpha in degrees(LE) xLE (from prev TE)
    % yLE (from prev TE)
    
    
    nacaCell=regexp(nacaStr,';','split');
    nacaCell=regexp(nacaCell,'_','split');
    
    nPtsPloop=nPtsPloop/str2double(nacaCell{1}{1});
    
    x=[linspace(0,1,round(nPtsPloop/2))];
    x=(0.5-cos(x*pi)/2);
    TEPos=[0 0];
    for ii=1:str2double(nacaCell{1}{1})
        refPos=TEPos+str2num([nacaCell{ii+1}{4},' ',nacaCell{ii+1}{5}]);
        [nacaCoord(ii).coord,TEPos]=GenerateNACACoordOrient...
            (x,nacaCell{ii+1}{1},nacaCell{ii+1}{2},nacaCell{ii+1}{3},...
            refPos,naca45t,naca4c,a4_closed,teps);
    end
    
    normL=sqrt(sum(TEPos.^2));
    rot=-asin(TEPos(2)/normL);
    for ii=1:str2double(nacaCell{1}{1})
        
        nacaCoord(ii).coord=nacaCoord(ii).coord/normL;
        rotMat=[cos(rot),-sin(rot);sin(rot),cos(rot)];
        nacaCoord(ii).coord=(rotMat*(nacaCoord(ii).coord)')';
    end
    
    
end

function [coord,TEPos]=GenerateNACACoordOrient(x,naca4Num,l,rot,refPos,tFunc,cFunc,a4,teps)
    
    [ctc,pct,tmax,refFlag]=ReadNacaString(naca4Num);
    
    tDist=tFunc(x,tmax,1,0,a4,teps);
    cDist=cFunc(x,ctc,pct,1,0)';
    
    unitCoord=RemoveIdenticalConsecutivePoints([[x',(tDist'+cDist)];flip([x',(-tDist'+cDist)])]);
    
    unitCoord=unitCoord*str2double(l)/10;
    [~,iTE]=max(unitCoord(:,1));
    rot=-str2double(rot)/180*pi;
    rotMat=[cos(rot),-sin(rot);sin(rot),cos(rot)];
    
    rotCoord=(rotMat*(unitCoord)')';
    
    coord=rotCoord+repmat(refPos,[size(rotCoord,1),1]);
    TEPos=coord(iTE,:);
end

function [airfoilDat,airfoil]=ReadAirfoilData(airfoilstr,airfoilDir)
    
    if nargin==1 || isempty(airfoilDir)
        airfoilDir='..\AeroLib\Smoothed\';
    end
    
    fid=fopen([airfoilDir,airfoilstr,'.dat'],'r');
    if fid==-1
        error('Airfoil Data could not be read. %s was not recognised as a valid file',[airfoilDir,airfoilstr,'.dat']);
    end
    for ii=1:3, str=fgetl(fid);end
    
    [airfoilDat,check]=fscanf(fid,'%f %f \n',[2,inf]);
    airfoilDat=airfoilDat';
    fclose(fid);
    
    [isCCW]=CCWLoop(airfoilDat);
    if ~isCCW
        airfoilDat=flip(airfoilDat,1);
    end
    iTE=find(airfoilDat(:,1)==max(airfoilDat(:,1)));
    if numel(iTE)==1
        airfoilDat=airfoilDat([iTE:end,1:iTE],:); % repeats the trailing edge at start and end.
    else
        [~,teT]=max(airfoilDat(iTE,2));
        [~,teB]=min(airfoilDat(iTE,2));
        if iTE(teT)==1
            airfoilDat=airfoilDat([iTE(teT):iTE(teB)],:);
        else
            airfoilDat=airfoilDat([iTE(teT):end,1:iTE(teB)],:);
        end
    end
    %analysisCoord(:,2)=analysisCoord(:,2)-analysisCoord(1,2); % 0 the y at the trailing edge.
    % TE is start and end of coord list, need LE
    %[dLE,iLE]=max(sum((analysisCoord-ones([size(analysisCoord,1),1])*analysisCoord(1,:)).^2,2));
    [~,iLE]=min(airfoilDat(:,1));
    if numel(iLE)==1
        airfoilDat=airfoilDat([1:iLE,iLE:end],:); % repeats the trailing edge at start and end.
        leT=1;
        leB=2;
        iLE=[iLE,iLE+1];
    else
        [~,leT]=max(airfoilDat(iLE,2));
        [~,leB]=min(airfoilDat(iLE,2));
        airfoilDat=airfoilDat([1:iLE(leT),iLE(leB):end],:);
    end
    
    
    upperLower=ones(size(airfoilDat(:,1)));
    teps=ones(size(airfoilDat(:,1)))*airfoilDat(iTE(teT),2);
    upperLower(iLE(leB):end)=-1;
    teps(iLE(leB):end)=airfoilDat(iTE(teB),2);
    airfoil.func.xMin=min(airfoilDat(:,1));
    airfoil.func.xMax=max(airfoilDat(:,1));
    airfoilDat(:,2)=airfoilDat(:,2)-teps.*(airfoilDat(:,1)-airfoil.func.xMin)...
        /(airfoil.func.xMax-airfoil.func.xMin);
    
    airfoil.coord=airfoilDat;
    airfoil.upperlower=upperLower;
    airfoil.func.upper=griddedInterpolant(flip(airfoilDat(1:iLE(leT),1)),...
        flip(airfoilDat(1:iLE(leT),2)),'linear','none');
    airfoil.func.lower=griddedInterpolant(airfoilDat(iLE(leB):end,1),...
        airfoilDat(iLE(leB):end,2),'linear','none');
end

function [constrVal]=OuterLimitInverse(gridrefined,paramoptim,airfoilstr)
    
    
    
    integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*...
        (tDistrib(1:end-1)+tDistrib(2:end))/2]);
    teps=5.48e-04/2/0.8;
    
    
    warning('Will only work with square or rectangular grids')
    varExtract={'axisRatio'};
    axisRatio=ExtractVariables(varExtract,paramoptim.parametrisation);
    
    [airfoilDat,airfoil]=ReadAirfoilData(airfoilstr,'');
    
    
    cellCentredGrid=CellCentreGridInformation(gridrefined);
    isActive=logical([cellCentredGrid(:).isactive]);
    %     actVerts=[cellCentredGrid(isActive).vertex];
    %     coord=vertcat(actVerts(:).coord);
    %     figure, hold on
    %     coord(:,2)=coord(:,2)*axisRatio;%;
    
    xMax=airfoil.func.xMax; %max(xPos);
    xMin=airfoil.func.xMin; %min(xPos);
    fillSub=zeros([1,sum(isActive)]);
    reqFrac=zeros([1,sum(isActive)]);
    actCellSub=find(isActive);
    for ii=1:numel(actCellSub)
        
        cellCoords=vertcat(cellCentredGrid(actCellSub(ii)).vertex(:).coord);
        cellCoords(:,2)=cellCoords(:,2)*axisRatio;
        [cornerCoord]=IdentifyCorners(cellCoords);
        posMin=min(cornerCoord);
        posMax=max(cornerCoord);
        
        %x=linspace(posMin(1),posMax(1),200);
        x=min(max(linspace(posMin(1),posMax(1),200),xMin),xMax);
        uppersurf=airfoil.func.upper(x)+(x-xMin)/(xMax-xMin)*teps;
        lowersurf=airfoil.func.lower(x)-(x-xMin)/(xMax-xMin)*teps;
        y=min(max(uppersurf,posMin(2)),posMax(2))-min(max(lowersurf,posMin(2)),posMax(2));
        %plot(x,cDistrib+y+posMin(2))
        vol=integr(x,y);
        fillSub(ii)=ii;
        reqFrac(ii)=vol(end)/cellCentredGrid(actCellSub(ii)).volume/axisRatio;
    end
    constrVal={fillSub,reqFrac};
end

function [constrVal]=MatchVoltoShape(gridrefined,paramoptim,shapePath)
    
    
    %     ((x(x>=(p*c))-xMin)/c)
    %
    %     [m*x(x<(cp*c))/p^2.*(2*p-x(x<(p*c))) ;...
    %         m*(1-x(x>=(p*c))/c)/(1-p)^2.*(1+x(x>=(p*c))-2*p)];
    
    integr=@(x,tDistrib) cumsum([0,(-x(1:end-1)+x(2:end)).*...
        (tDistrib(1:end-1)+tDistrib(2:end))/2]);
    
    warning('Will only work with square or rectangular grids')
    varExtract={'axisRatio'};
    axisRatio=ExtractVariables(varExtract,paramoptim.parametrisation);
    [uppersurf,lowersurf]=ReadShapeIn(shapePath);
    
    
    cellCentredGrid=CellCentreGridInformation(gridrefined);
    isActive=logical([cellCentredGrid(:).isactive]);
    actVerts=[cellCentredGrid(isActive).vertex];
    coord=vertcat(actVerts(:).coord);
    %figure, hold on
    coord(:,2)=coord(:,2)*axisRatio;%;
    
    xPos=RemoveIdenticalEntries(coord(:,1));
    xMax=max(xPos);
    xMin=min(xPos);
    xDist=(xMax-xMin);
    fillSub=zeros([1,sum(isActive)]);
    reqFrac=zeros([1,sum(isActive)]);
    actCellSub=find(isActive);
    for ii=1:numel(actCellSub)
        
        cellCoords=vertcat(cellCentredGrid(actCellSub(ii)).vertex(:).coord);
        cellCoords(:,2)=cellCoords(:,2)*axisRatio;
        [cornerCoord]=IdentifyCorners(cellCoords);
        posMin=min(cornerCoord);
        posMax=max(cornerCoord);
        
        x=linspace(posMin(1),posMax(1),200);
        
        yup=interp1(uppersurf(:,1)*xDist+xMin,uppersurf(:,2)*xDist,x);
        ylow=interp1(lowersurf(:,1)*xDist+xMin,lowersurf(:,2)*xDist,x);
        y=min(max(yup,posMin(2)),posMax(2))-min(max(ylow,posMin(2)),posMax(2));
        %plot(x,cDistrib+y+posMin(2))
        vol=integr(x,y);
        fillSub(ii)=ii;
        reqFrac(ii)=vol(end)/cellCentredGrid(actCellSub(ii)).volume/axisRatio;
    end
    constrVal={fillSub,reqFrac};
end

function [uppersurf,lowersurf]=ReadShapeIn(shapepath)
    
    filename=dir(shapepath);
    filename=filename(1).name;
    extPos=regexp(filename,'\.');
    ext=filename(extPos+1:end);
    
    switch ext
        case 'mat'
            load(shapepath)
            if ~exist('uppersurf','var');error('Data file loaded did not have the upper surface');end
            if ~exist('lowersurf','var');error('Data file loaded did not have the lower surface');end
            
        case 'dat'
            error('not coded')
            
        otherwise
            error('Unknown type')
            
    end
    
end

function [ctc,pct,ttc,refFlag]=ReadNacaString(nacaStr)
    refFlag=0;
    if numel(nacaStr)==4
        ctc=str2num(nacaStr(1))/100;
        pct=str2num(nacaStr(2))/10;
        ttc=str2num(nacaStr(3:4))/100;
        
        
    elseif numel(nacaStr)==5
        ctc=str2num(nacaStr(1))/100;
        pct=str2num(nacaStr(2))/10;
        refFlag=str2num(nacaStr(3));
        ttc=str2num(nacaStr(4:5))/100;
        %error('Five digits not implemented - need for tabulated data')
    else
        error('Invalid string of length %i',numel(nacaStr))
    end
end

function [cornersCoord]=IdentifyCorners(coord)
    
    minXPos=find(coord(:,1)==min(coord(:,1)));
    maxXPos=find(coord(:,1)==max(coord(:,1)));
    minYPos=find(coord(:,2)==min(coord(:,2)));
    maxYPos=find(coord(:,2)==max(coord(:,2)));
    
    [~,loleCorn]=min(coord(minXPos,2));
    loleCorn=minXPos(loleCorn);
    [~,hileCorn]=max(coord(minXPos,2));
    hileCorn=minXPos(hileCorn);
    [~,loriCorn]=min(coord(maxXPos,2));
    loriCorn=maxXPos(loriCorn);
    [~,hiriCorn]=max(coord(maxXPos,2));
    hiriCorn=maxXPos(hiriCorn);
    
    
    [~,leloCorn]=min(coord(minYPos,1));
    leloCorn=minYPos(leloCorn);
    [~,riloCorn]=max(coord(minYPos,1));
    riloCorn=minYPos(riloCorn);
    [~,lehiCorn]=min(coord(maxYPos,1));
    lehiCorn=maxYPos(lehiCorn);
    [~,rihiCorn]=max(coord(maxYPos,1));
    rihiCorn=maxYPos(rihiCorn);
    
    tests=(loleCorn==leloCorn) && (loriCorn==riloCorn) &&(hileCorn==lehiCorn) && (hiriCorn==rihiCorn);
    if tests==false
        warning('Unexpected Corners returning only one set')
    end
    cornersCoord=coord([loleCorn,loriCorn,hileCorn,hiriCorn],:);
end

function [desVarList]=ExtractActiveVariable(nFill,notDesInd,inactiveVar)
    
    desVarList=1:nFill;
    desVarList([notDesInd,inactiveVar])=[];
    
end

function [inactiveVar]=SelectInactiveVariables(newFill,varActive,derivtenscalc,cutoff)
    
    switch varActive
        case 'all'
            inactiveVar=[];
        case 'border'
            [inactiveVar]=InactiveVariables_border(newFill);
            
        case 'wideborder'
            [inactiveVar,activeVar]=InactiveVariables_wideborder(newFill,derivtenscalc,cutoff);
            
        case 'snaksensiv'
            
            [inactiveVar]=InactiveVariables_border(newFill);
        otherwise
            error('unrecognised variable activation criterion')
    end
    
    
end

function [inactiveVar,activeVar]=InactiveVariables_border(newFill)
    
    is0=newFill==0;
    is1=newFill==1;
    
    inactiveVar=find(is0);
    inactiveVar=[inactiveVar,find(is1)];
    
    isBord=~is0 & ~is1;
    activeVar=find(isBord);
    
    
end

function [inactiveVar,activeVar]=InactiveVariables_wideborder(newFill,derivtenscalc,cutoff)
    
    [inactiveVar,activeVar]=InactiveVariables_border(newFill);
    
    indexDeriv=[derivtenscalc(:).index];
    actVarSub=FindObjNum([],activeVar,indexDeriv);
    actFill=newFill(activeVar);
    cutActVar=actFill>cutoff;
    activeVar=unique([activeVar,[derivtenscalc(actVarSub(cutActVar)).neighbours]]);
    inactiveVar=1:length(newFill);
    inactiveVar(activeVar)=[];
end

function [isGradient]=CheckIfGradient(optimMethod)
    
    switch optimMethod
        
        case 'DE'
            isGradient=false;
        case 'DEtan'
            isGradient=false;
        case 'conjgrad'
            isGradient=true;
        case 'conjgradls'
            isGradient=true;
        otherwise
            isGradient=false;
            warning('Optimisation method is not known as gradient based or otherwise, no gradient is assumed')
            
    end
end

function [isAnalytical]=CheckIfAnalytical(objectiveName)
    
    switch objectiveName
        
        case 'Rosenbrock'
            isAnalytical=true;
        otherwise
            isAnalytical=false;
            warning('Objective is not specified assumed non analytical function')
            
    end
end

function [isSnakeSensitivity]=CheckSnakeSensitivityAlgorithm(paramoptim)
    
    varExtract={'varActive','optimMethod'};
    [varActive,optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    isSnakeSensitivity=false;
    [isGradient]=CheckIfGradient(optimMethod);
    
    if isGradient
        isSnakeSensitivity=strcmp(varActive,'snaksensiv');
    end
end

function [newRootFill,optionOut]=OverflowHandling(paramoptim,newRootFill,extraarguments)
    % Function which handles steps overflowing the fill constraint
    
    varExtract={'desVarRange','varOverflow','nDesVar'};
    [desVarRange,varOverflow,nDesVar]=ExtractVariables(varExtract,paramoptim);
    optionOut=[];
    switch varOverflow
        case 'truncate'
            minD=min(desVarRange);
            maxD=max(desVarRange);
            newRootFill(newRootFill<minD)=minD;
            newRootFill(newRootFill>maxD)=maxD;
            
        case 'spill'
            newFill=zeros(size(newRootFill));
            for ii=1:size(newRootFill,1)
                [newFill(ii,1:nDesVar)]=SpillOverflow(paramoptim,newRootFill(ii,1:nDesVar));
            end
            newRootFill(:,1:nDesVar)=newFill(:,1:nDesVar);
            minD=min(desVarRange);
            maxD=max(desVarRange);
            newRootFill(newRootFill<minD)=minD;
            newRootFill(newRootFill>maxD)=maxD;
        case 'vertexflow'
            
            if nargin>2
                try
                    [newRootFill(:,1:nDesVar),optionOut]=RunVertexOverflow(paramoptim,newRootFill(:,1:nDesVar),extraarguments);
                    disp('Vertex Flow succesful')
                catch ME
                    disp(ME.getReport)
                    warning('Vertex Flow Failed')
                end
            end
            
            newFill=zeros(size(newRootFill));
            for ii=1:size(newRootFill,1)
                [newFill(ii,1:nDesVar)]=SpillOverflow(paramoptim,newRootFill(ii,1:nDesVar));
            end
            newRootFill(:,1:nDesVar)=newFill(:,1:nDesVar);
            minD=min(desVarRange);
            maxD=max(desVarRange);
            newRootFill(newRootFill<minD)=minD;
            newRootFill(newRootFill>maxD)=maxD;
        otherwise
            error('No valid design variable overflow mechanism')
            
    end
    
    
    
end

function [newFill]=SpillOverflow(paramoptim,newRootFill)
    
    varExtract={'desVarRange','desvarconnec'};
    [desVarRange,desvarconnec]=ExtractVariables(varExtract,paramoptim);
    
    normRootFill=(newRootFill-min(desVarRange))/(max(desVarRange)-min(desVarRange));
    
    overVar=find(normRootFill>1);
    underVar=find(normRootFill<0);
    exitFlag=false;
    kk=0;
    while (~isempty(overVar) || ~isempty(underVar)) && ~exitFlag
        
        
        [normRootFill]=SpillOverflowVarHandling(normRootFill,desvarconnec,overVar);
        
        normRootFill=1-normRootFill;
        [normRootFill]=SpillOverflowVarHandling(normRootFill,desvarconnec,underVar);
        normRootFill=1-normRootFill;
        
        
        
        overVar=find(normRootFill>1);
        underVar=find(normRootFill<0);
        meanFill=mean(normRootFill);
        stdFill=std(normRootFill);
        kk=kk+1;
        if ((meanFill>=1 || meanFill<=0)  && (abs(stdFill)<1e-5)) || kk>100
            exitFlag=true;
        end
        
    end
    
    newFill=normRootFill*(max(desVarRange)-min(desVarRange))+min(desVarRange);
end

function [newRootFill]=SpillOverflowVarHandling(newRootFill,desvarconnec,flowVar)
    
    desVarIndList=[desvarconnec(:).index];
    
    overFlowMat=zeros([length(flowVar),length(newRootFill)]);
    
    for ii=1:length(flowVar)
        
        currSub=FindObjNum([],flowVar(ii),desVarIndList);
        
        neighInd=desvarconnec(currSub).neighbours;
        cornInd=desvarconnec(currSub).corners;
        
        neighSub=FindObjNum([],neighInd,desVarIndList);
        cornSub=FindObjNum([],cornInd,desVarIndList);
        
        currVol=newRootFill(currSub);
        neighVol=newRootFill(neighSub);
        cornVol=newRootFill(cornSub);
        
        neighEmpt=neighSub(neighVol<=0);
        cornEmpt=cornSub(cornVol<=0);
        neighGreyCell={[]};
        for jj=1:length(cornEmpt)
            neighGreyCell{jj}=FindObjNum([],desvarconnec(cornEmpt(jj)).neighbours,neighInd)';
        end
        neighAll=[neighGreyCell{:}];
        neighSubAll=neighSub(RemoveIdenticalEntries(neighAll(neighAll~=0)));
        neighSubAll=RemoveIdenticalEntries([neighEmpt;neighSubAll]);
        neighSubAll=neighSubAll(newRootFill(neighSubAll)<1);
        
        nCorn=numel(cornEmpt);
        nNeigh=sum(1-newRootFill(neighSubAll));
        isNormProb=true;
        
        if nCorn>0
            baseRate=(-nNeigh+sqrt(nNeigh^2+4*nCorn))/(2*nCorn);
            
        elseif nNeigh>0
            baseRate=1/nNeigh;
            
        elseif numel(neighSub(newRootFill(neighSub)<1))>0 ...
                && numel(cornSub(newRootFill(cornSub)<1))>0
            
            neighSubAll=neighSub(newRootFill(neighSub)<1);
            cornEmpt=cornSub(newRootFill(cornSub)<1);
            nCorn=numel(cornEmpt);
            nNeigh=sum(1-newRootFill(neighSubAll));
            
            baseRate=(-nNeigh+sqrt(nNeigh^2+4*nCorn))/(2*nCorn);
            
        elseif numel(neighSub(newRootFill(neighSub)<1))>0
            neighSubAll=neighSub(newRootFill(neighSub)<1);
            nNeigh=sum(1-newRootFill(neighSubAll));
            baseRate=1/nNeigh;
            
        elseif numel(cornSub(newRootFill(cornSub)<1))>0
            neighSubAll=cornSub(newRootFill(cornSub)<1);
            nNeigh=sum(1-newRootFill(neighSubAll));
            baseRate=1/nNeigh;
            
        elseif numel(neighSub(newRootFill(neighSub)<currVol))>0 || ...
                numel(cornSub(newRootFill(cornSub)<currVol))>0
            neighSubAll=[neighSub(newRootFill(neighSub)<currVol);cornSub(newRootFill(cornSub)<currVol)];
            nNeigh=sum(currVol-newRootFill(neighSubAll));
            baseRate=1/nNeigh;
            
            isNormProb=false;
        else
            neighSubAll=[];
            nNeigh=sum(-newRootFill(neighSubAll));
            
            baseRate=0;
            
            isNormProb=false;
        end
        
        if isNormProb
            overFlowVol=zeros(size(newRootFill));
            overFlowVol(neighSubAll)=(1-newRootFill(neighSubAll))*baseRate;
            overFlowVol(cornEmpt)=baseRate^2;
            overFlowVol(currSub)=-1;
            overFlowMat(ii,:)=overFlowVol*(currVol-1);
            %newRootFill=newRootFill+overFlowVol*(currVol-1);
        else
            overFlowVol=zeros(size(newRootFill));
            overFlowVol(neighSubAll)=baseRate*(currVol-newRootFill(neighSubAll));
            overFlowVol(cornEmpt)=baseRate^2;
            overFlowVol(currSub)=-sum(overFlowVol);
            overFlowMat(ii,:)=overFlowVol*(currVol-max([min(newRootFill(neighSubAll)),1]))*2/3;
            %newRootFill=newRootFill+overFlowVol*(currVol-max([min(newRootFill(neighSubAll)),1]))*2/3;
        end
        
    end
    
    newRootFill=newRootFill+sum(overFlowMat,1);
    
end
%% Vertex flow functions

function [newFill,fillflow]=RunVertexOverflow(paramoptim,newFill,pathSnak)
    if ~isstruct(pathSnak)
        [gridBase,gridRefined,gridConnec,cellCentredGrid,snaxel]=ExtractSnakInfo_vertflow(pathSnak);
        [fillflow]=VertexOverFlowDefine(paramoptim,gridBase,gridRefined,gridConnec,cellCentredGrid,snaxel);
    else
        fillflow=pathSnak;
    end
    for ii=1:size(newFill,1)
        [newFill(ii,:)]=VertexOverFlowExecute(fillflow,newFill(ii,:));
        
    end
end

function [gridBase,gridRefined,gridConnec,cellCentredGrid,snaxel]=ExtractSnakInfo_vertflow(pathSnak)
    
    [returnPath,returnName]=FindDir(pathSnak,'restart',0);
    
    load(returnPath{1},'restartsnak')
    cellCentredGrid=restartsnak.cellCentredGrid;
    snaxel=restartsnak.snaxel;
    
    gridPath=[regexprep(returnPath{1},returnName{1},''),filesep,'..',filesep,...
        '..',filesep,'iteration_0',filesep,'profile_0',filesep];
    
    [returnPath,returnName]=FindDir(gridPath,'restart',0);
    load(returnPath{1},'grid')
    gridBase=grid.base;
    gridRefined=grid.refined;
    gridConnec=grid.connec;
end

function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
    %     if iscell(rootDir)
    %         subDir=dir(rootDir{1});
    %         subDir(1:2)=[];
    %         for ii=2:numel(rootDir)
    %             partsubDir=dir(rootDir{ii});
    %             partsubDir(1:2)=[];
    %             subDir=[subDir,partsubDir];
    %         end
    %     else
    subDir=dir(rootDir);
    subDir(1:2)=[];
    %     end
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    
    
    if isempty(returnSub)
        fprintf('FindDir Could not find requested item %s in:\n%s \n',strDir,rootDir)
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
    
    
end

function [fillflow]=VertexOverFlowDefine(paramoptim,gridBase,gridRefined,gridConnec,cellCentredGrid,snaxel)
    
    snakGridRef=ExtractVariables({'refineGrid'},paramoptim.parametrisation);
    
    try
        varExtract={'spillCutOff'};
        [spillCutOff]=ExtractVariables(varExtract,paramoptim);
    catch ME
        disp(ME.getReport)
        spillCutOff=2e-2;
    end
    
    [snaxCorn,snaxBord]=ExtractVerticesForFlow(gridBase,gridRefined,gridConnec,snaxel,snakGridRef);
    snaxFlow=[snaxCorn;snaxBord];
    
    vertAct=[cellCentredGrid([cellCentredGrid(:).fill]==0 | [cellCentredGrid(:).fill]==1).vertex];
    vertAct=unique([vertAct(:).index]);
    
    snaxFlow=snaxFlow(find(FindObjNum([],snaxFlow(:,1),vertAct)~=0),:);
    
    % Fill dist provides an indication of the fill distance rather than the
    % actual
    if ~isempty(snaxFlow)
        fillDist=unique(snaxFlow(:,1));
        for ii=1:numel(fillDist)
            fillDist(ii,2)=prod(snaxFlow(find(snaxFlow(:,1)==fillDist(ii,1)),2));
        end
        fillDist=fillDist(find(fillDist(:,2)<spillCutOff),:);
        
        % need to transform fill dist into a relationship between cells (ie
        % filling cells and emptying ones)
        
        vertflowinfo=BuildVertexOverflowStruct(fillDist,gridRefined,cellCentredGrid);
        fillflow=BuildFillFlow(vertflowinfo,cellCentredGrid,gridRefined,gridBase,gridConnec);
    else
        fillflow=struct([]);
    end
end

function [vertflowinfo]=BuildVertexOverflowStruct(fillDist,gridRefined,cellCentredGrid)
    
    
    edgeVertInd=[gridRefined.edge(:).vertexindex];
    vertInd=[gridRefined.vertex(:).index];
    cellInd=[cellCentredGrid(:).index];
    vertflowinfo=repmat(struct('vertex',[],'cells',repmat(struct('index',[],'edges',[],...
        'issource',[]),[1 0]),'filldist',[]),[1,size(fillDist,1)]);
    
    for ii=1:size(fillDist,1)
        
        vertflowinfo(ii).vertex=fillDist(ii,1);
        vertflowinfo(ii).filldist=fillDist(ii,2);
        currEdges=ceil(FindObjNum([],fillDist(ii,1),edgeVertInd)/2);
        
        currVertind=[gridRefined.edge(currEdges).vertexindex];
        currVertind=currVertind(currVertind~=fillDist(ii,1));
        currCoords=vertcat(gridRefined.vertex(FindObjNum([],[fillDist(ii,1),currVertind],vertInd)).coord);
        currVecs=currCoords(2:end,:)-(ones([size(currCoords,1)-1,1]))*currCoords(1,:);
        
        [~,iVec]=sort(ExtractAngle360(currVecs(1,:),currVecs));
        currVecs=currVecs(iVec,:);
        currEdges=currEdges(iVec);
        currEdgesInd=[gridRefined.edge(currEdges).index];
        
        [currcells,~,cells2edges]=unique([gridRefined.edge(currEdges).cellindex]);
        
        for jj=numel(currcells):-1:1;
            
            vertflowinfo(ii).cells(jj).index=currcells(jj);
            c2e=ceil(find(cells2edges==jj)/2);
            vertflowinfo(ii).cells(jj).edges=currEdgesInd(c2e);
            
            vertflowinfo(ii).cells(jj).angleratio=abs(ExtractAngle360(...
                currVecs(c2e(1),:),currVecs(c2e(2),:))/2/pi);
            if abs(c2e(2)-c2e(1))>1
                vertflowinfo(ii).cells(jj).angleratio=1-vertflowinfo(ii).cells(jj).angleratio;
            end
            
            currSub=FindObjNum([],currcells(jj),cellInd);
            if currSub~=0
                vertflowinfo(ii).cells(jj).issource=~(cellCentredGrid(currSub).fill==0 ...
                    || cellCentredGrid(currSub).fill==1);
            end
        end
        
    end
    
end

function [fillflow]=BuildFillFlow(vertflowinfo,cellCentredGrid,gridRefined,gridBase,gridConnec)
    
    oldCellInd=ReverseStructInfo(gridConnec.cell,'oldCellInd','newCellInd');
    newCellInd=[gridConnec.cell(:).newCellInd];
    cellInd=[cellCentredGrid(:).index];
    activeInd=[cellCentredGrid(logical([cellCentredGrid(:).isactive])).index];
    % Normalise and pass the data to the large cell
    for ii=1:numel(vertflowinfo)
        
        currOld=oldCellInd(FindObjNum([],[vertflowinfo(ii).cells(:).index],newCellInd));
        rmRow=false(size(currOld));
        for jj=1:numel(vertflowinfo(ii).cells)
            vertflowinfo(ii).cells(jj).index=currOld(jj);
            firstMatch=min(find(currOld==currOld(jj)));
            rmRow(jj)=firstMatch<jj;
            if rmRow(jj)
                vertflowinfo(ii).cells(firstMatch).angleratio=...
                    vertflowinfo(ii).cells(firstMatch).angleratio...
                    +vertflowinfo(ii).cells(jj).angleratio;
                
                tempEdge=sort([vertflowinfo(ii).cells(firstMatch).edges,...
                    vertflowinfo(ii).cells(jj).edges]);
                iskeep=~((tempEdge-tempEdge([2:end,1]))==0 | (tempEdge-tempEdge([end,1:end-1]))==0);
                vertflowinfo(ii).cells(firstMatch).edges=tempEdge(iskeep);
                
            end
        end
        vertflowinfo(ii).cells(rmRow)=[];
    end
    for ii=1:numel(vertflowinfo)
        angleRat=[vertflowinfo(ii).cells(:).angleratio];
        isSource=[vertflowinfo(ii).cells(:).issource];
        angleRat(isSource)=angleRat(isSource)/sum(angleRat(isSource));
        angleRat(~isSource)=angleRat(~isSource)/sum(angleRat(~isSource));
        
        
        for jj=1:numel(vertflowinfo(ii).cells)
            vertflowinfo(ii).cells(jj).angleratio=angleRat(jj);
        end
    end
    
    %
    fillflow=repmat(struct('cellsource',[],'celldest',[],'fillcutoff',[],'floodratio',[],'dir',[]),[1 0]);
    kk=0;
    for ii=1:numel(vertflowinfo)
        
        for jj=find([vertflowinfo(ii).cells(:).issource])
            for ll=find(~[vertflowinfo(ii).cells(:).issource])
                kk=kk+1;
                fillflow(kk).cellsource=vertflowinfo(ii).cells(jj).index;
                fillflow(kk).celldest=vertflowinfo(ii).cells(ll).index;
                sourceSub=FindObjNum([],fillflow(kk).cellsource,cellInd);
                destSub=FindObjNum([],fillflow(kk).celldest,cellInd);
                dir=((cellCentredGrid(destSub).fill<0.5)-0.5)*2;
                fillflow(kk).fillcutoff=cellCentredGrid(sourceSub).fill...
                    +dir*vertflowinfo(ii).filldist;
                fillflow(kk).dir=dir;
                
                fillflow(kk).floodratio=CalculateFloodRatio(cellCentredGrid(sourceSub)...
                    ,cellCentredGrid(destSub),vertflowinfo(ii).cells([jj,ll]));
                fillPos=FindObjNum([],[fillflow(kk).cellsource,fillflow(kk).celldest],activeInd);
                if any(fillPos==0)
                    fillflow(kk)=[];
                    kk=kk-1;
                else
                    fillflow(kk).cellsource=fillPos(1);
                    fillflow(kk).celldest=fillPos(2);
                end
            end
        end
        
        
    end
    
    delFillFlow=false(size(fillflow));
    fillSource=[fillflow(:).cellsource];
    fillDest=[fillflow(:).celldest];
    for ii=1:numel(fillflow)
        if ~delFillFlow(ii)
            
            isSame=find(fillSource(ii)==fillSource & fillDest(ii)==fillDest);
            if numel(isSame)>1
                if fillflow(ii).dir
                    fillCutOff=min([fillflow(isSame).fillcutoff]);
                    floodRatio=max([fillflow(isSame).floodratio]);
                else
                    fillCutOff=max([fillflow(isSame).fillcutoff]);
                    floodRatio=max([fillflow(isSame).floodratio]);
                end
                fillflow(isSame(1)).fillcutoff=fillCutOff;
                fillflow(isSame(1)).floodratio=floodRatio;
                delFillFlow(isSame(2:end))=true;
            end
            
        end
    end
    fillflow(delFillFlow)=[];
end

function [cellRatio]=CalculateFloodRatio(cellSource,cellDest,cellflowinfo)
    % calculates the bleed ratio between the source and destination cell
    % based on wether it is a corner bleed and
    cellRatio=1/4*cellSource.volume/cellDest.volume;
    
end

function [snaxCorn,snaxBord]=ExtractVerticesForFlow(gridBase,gridRefined,...
        gridConnec,snaxel,snakGridRef)
    % Finds vertices where the snake is close to the edge of a design
    % variable and returns the concerned vertices along with the normalised
    % length away from it.
    
    
    vertInd=[gridBase.vertex(:).index];
    if numel(snakGridRef)>1
        
        edgeInd=[gridRefined.edge(:).index];
        edgeOrient=[gridRefined.edge(:).orientation];
        edgeOrient=edgeOrient(FindObjNum([],[snaxel(:).edge],edgeInd))+1;
        snaxD=[[[snaxel(:).tovertex];(1-[snaxel(:).d])./snakGridRef(edgeOrient)],...
            [[snaxel(:).fromvertex];[snaxel(:).d]./snakGridRef(edgeOrient)]]';
    else
        
        snaxD=[[[snaxel(:).tovertex];(1-[snaxel(:).d])/snakGridRef(1)],...
            [[snaxel(:).fromvertex];[snaxel(:).d]/snakGridRef(1)]]';
    end
    snaxVertSub=FindObjNum([],snaxD(:,1),vertInd);
    snaxCorn=snaxD(find(snaxVertSub~=0),:);
    snaxCornVertS=snaxVertSub(snaxVertSub~=0);
    %
    
    oldCellDat=[gridConnec.edge(:).newcell];
    [oldIndsNewOrd]=[0,ReverseStructInfo(gridConnec.cell,'oldCellInd','newCellInd')];
    oldCellDat=oldIndsNewOrd(FindObjNum([],oldCellDat,[0,gridConnec.cell(:).newCellInd]));
    borderedges=[gridConnec.edge(oldCellDat(1:2:end)~=oldCellDat(2:2:end)).newedge];
    borderVertex=unique([gridRefined.edge(FindObjNum([],borderedges,...
        [gridRefined.edge(:).index])).vertexindex]);
    borderVertex(FindObjNum([],vertInd,borderVertex))=[];
    snaxVertSub=FindObjNum([],snaxD(:,1),borderVertex);
    snaxBord=snaxD(find(snaxVertSub~=0),:);
    snaxBordVertS=snaxVertSub(snaxVertSub~=0);
    
end

function [newFill]=VertexOverFlowExecute(fillflow,newFill)
    
    for ii=1:numel(fillflow)
        bleedOff=newFill(fillflow(ii).cellsource)-fillflow(ii).fillcutoff;
        if sign(bleedOff)~=(fillflow(ii).dir*2-1)
            bleedOff=0;
        end
        bleedOff=bleedOff*(fillflow(ii).floodratio);
        newFill(fillflow(ii).celldest)=newFill(fillflow(ii).celldest)+bleedOff;
    end
    
    
end

%%

function population=ApplySymmetry(paramoptim,population)
    
    varExtract={'symDesVarList'};
    [symDesVarList]=ExtractVariables(varExtract,paramoptim);
    
    for ii=1:length(population)
        population(ii).fill(symDesVarList(2,:))=...
            population(ii).fill(symDesVarList(1,:));
    end
    
    
end

function [res]=TestDeviationPopStruct(population,cellTest,dirOut)
    
    for ii=1:numel(cellTest)
        analysisVal=[];
        for jj=1:numel(population)
            analysisVal=[analysisVal,population(jj).additional.(cellTest{ii})];
        end
        h(ii)=figure('Name',cellTest{ii});
        subplot(2,2,1)
        normplot(analysisVal)
        res.(cellTest{ii}).mean=mean(analysisVal);
        res.(cellTest{ii}).median=median(analysisVal);
        res.(cellTest{ii}).std=std(analysisVal);
        res.(cellTest{ii}).max=max(analysisVal);
        res.(cellTest{ii}).min=min(analysisVal);
        res.(cellTest{ii}).decile=quantile(analysisVal,19);
        res.(cellTest{ii}).ratiosuccess=sum(analysisVal<(8e-4))/numel(analysisVal);
        subplot(2,2,3)
        semilogy(res.(cellTest{ii}).decile)
        grid on
        xlabel('decile')
        ylabel('Val')
        analysisVal=log10(analysisVal);
        subplot(2,2,2)
        normplot(analysisVal)
        res.(cellTest{ii}).logmean=10^mean(analysisVal);
        res.(cellTest{ii}).logmedian=10^median(analysisVal);
        res.(cellTest{ii}).logstd=std(analysisVal);
        res.(cellTest{ii}).logmax=10^max(analysisVal);
        res.(cellTest{ii}).logmin=10^min(analysisVal);
        res.(cellTest{ii}).logdecile=10.^quantile(analysisVal,19);
        subplot(2,2,4)
        plot(quantile(analysisVal,19))
        grid on
        xlabel('decile')
        ylabel('log10(Val)')
    end
    
    if nargin>2
        for ii=1:numel(h)
            figName=[dirOut,filesep,'Figure_',int2str(ii),'_',h(ii).Name,'.fig'];
            hgsave(h(ii),figName);
        end
        fileOut=fopen([dirOut,filesep,'TabResSummary.csv'],'w');
        fieldsRes=fieldnames(res.(cellTest{1}));
        
        fprintf(fileOut,' ');
        for ii=1:numel(fieldsRes)
            fprintf(fileOut,', %s',fieldsRes{ii});
        end
        
        for jj=1:numel(cellTest)
            fprintf(fileOut,'\n');
            fprintf(fileOut,'%s',cellTest{jj});
            for ii=1:numel(fieldsRes)
                fprintf(fileOut,',%s',...
                    num2str(res.(cellTest{jj}).(fieldsRes{ii}),'%20.8e '));
            end
        end
        fclose(fileOut);
    end
end

function [paramoptim,paramsnake]=HandleNonFillVar(member,paramoptim,paramsnake)
    
    varExtract={'nonFillVar','numNonFillVar'};
    
    [nonFillVar,numNonFillVar]=ExtractVariables(varExtract,paramoptim);
    
    if nargin<3
        paramsnake=paramoptim.parametrisation;
    end
    
    for ii=1:numel(nonFillVar)
        vS=1+sum(numNonFillVar(1:ii-1));
        vE=sum(numNonFillVar(1:ii));
        switch nonFillVar{ii}
            case 'axisratio'
                paramoptim.parametrisation=SetVariables({'axisRatio'},...
                    {member.nonfillvar(vS:vE)},paramoptim.parametrisation);
                
                paramsnake=SetVariables({'axisRatio'},...
                    {member.nonfillvar(vS:vE)},paramsnake);
            case 'alpha'
                paramoptim=SetVariables({'nAlpha'},...
                    {member.nonfillvar(vS:vE)},paramoptim);
            case 'mach'
                paramoptim=SetVariables({'nMach'},...
                    {member.nonfillvar(vS:vE)},paramoptim);
                
            otherwise
                warning('Unknown Design Variable specified')
        end
        
        
    end
    
    
end

function [population]=RescaleDesVarNoFill(scaleDir,paramoptim,population)
    
    varExtract={'desVarRange','desVarRangeNoFill','nonFillVar','numNonFillVar','nDesVar'};
    
    [desVarRange,desVarRangeNoFill,nonFillVar,numNonFillVar,nDesVar]=ExtractVariables(varExtract,paramoptim);
    
    
    switch scaleDir
        case 'tofill'
            % scale from the fill range to the actual range
            for ii=1:numel(population)
                kk=0;
                for jj=1:numel(nonFillVar)
                    desCurr=population(ii).nonfillvar(kk+1:kk+numNonFillVar(jj));
                    population(ii).nonfillvar(kk+1:kk+numNonFillVar(jj))=...
                        ((desCurr-min(desVarRangeNoFill{jj}))/...
                        (max(desVarRangeNoFill{jj})-min(desVarRangeNoFill{jj})))...
                        *(max(desVarRange)-min(desVarRange))+min(desVarRange);
                    
                    
                    isAct=(population(ii).optimdat.var>nDesVar+kk) & ...
                        (population(ii).optimdat.var<nDesVar+kk+numNonFillVar(jj));
                    kk=kk+numNonFillVar(jj);
                end
                if isstruct(population(ii).optimdat)
                    kk=0;
                    for jj=1:numel(nonFillVar)
                        isAct=(population(ii).optimdat.var>nDesVar+kk) & ...
                            (population(ii).optimdat.var<=nDesVar+kk+numNonFillVar(jj));
                        datCurr=population(ii).optimdat.value(isAct);
                        population(ii).optimdat.value(isAct)=...
                            ((datCurr)/...
                            (max(desVarRangeNoFill{jj})-min(desVarRangeNoFill{jj})))...
                            *(max(desVarRange)-min(desVarRange));
                        
                        kk=kk+numNonFillVar(jj);
                    end
                end
            end
        case 'tovar'
            % scale from the fill range to the actual range
            for ii=1:numel(population)
                kk=0;
                for jj=1:numel(nonFillVar)
                    desCurr=population(ii).nonfillvar(kk+1:kk+numNonFillVar(jj));
                    population(ii).nonfillvar(kk+1:kk+numNonFillVar(jj))=...
                        ((desCurr-min(desVarRange))/...
                        (max(desVarRange)-min(desVarRange)))...
                        *(max(desVarRangeNoFill{jj})-min(desVarRangeNoFill{jj}))...
                        +min(desVarRangeNoFill{jj});
                    kk=kk+numNonFillVar(jj);
                end
                if isstruct(population(ii).optimdat)
                    kk=0;
                    for jj=1:numel(nonFillVar)
                        isAct=(population(ii).optimdat.var>nDesVar+kk) & ...
                            (population(ii).optimdat.var<=nDesVar+kk+numNonFillVar(jj));
                        datCurr=population(ii).optimdat.value(isAct);
                        population(ii).optimdat.value(isAct)=...
                            ((datCurr)/(max(desVarRange)-min(desVarRange)))...
                            *(max(desVarRangeNoFill{jj})-min(desVarRangeNoFill{jj}));
                        
                        kk=kk+numNonFillVar(jj);
                    end
                end
                
            end
            
    end
    
end

%% Flow

function [loop]=ConstantArea_Busemann(xMin,xMax,A,M)
    
    h=2*A/2/(xMax-xMin);
    Dx=xMax-xMin;
    tanDel=2*h/Dx;
    points=[xMin,0;xMax,0;xMin+(xMax-xMin)/2,h];
    
    [np]=GoldenSection_maxTanDel(0,pi/2,1e-6,M);
    [B]=GoldenSection_FindB(0,np,1e-6,M,tanDel);
    
    H=tan(B)*Dx/2;
    offset=H+h;
    points2=points;
    points2(:,2)=-points2(:,2)+offset/2;
    
    points(:,2)=points(:,2)-offset/2;
    points2=points2(end:-1:1,:);
    loop(1).subdivision=[points];
    loop(1).isccw=true;
    loop(2).subdivision=[points2];
    loop(2).isccw=true;
    loop(1).subdivision(end+1,:)= loop(1).subdivision(1,:);
    loop(2).subdivision(end+1,:)= loop(2).subdivision(1,:);
end

function [tanDel]=CalcTanDel(M,B)
    
    tanDel=2/tan(B)*(M^2*(sin(B)^2)-1)/(M^2*(1.4+cos(2*B))+2);
    
end

function [np]=GoldenSection_maxTanDel(lb,hb,tol,M)
    
    gr=(sqrt(5)-1)/2;
    
    c=hb-gr * (hb-lb);
    d=lb+gr * (hb-lb);
    
    while abs(c-d)>tol
        
        fc=-CalcTanDel(M,c);
        fd=-CalcTanDel(M,d);
        
        if fc<fd
            hb=d;
            
        else
            lb=c;
        end
        
        c=hb-gr * (hb-lb);
        d=lb+gr * (hb-lb);
    end
    
    np=(hb+lb)/2;
    
end

function [np]=GoldenSection_func(lb,hb,tol,extraIn,func)
    % finds a maximum
    
    gr=(sqrt(5)-1)/2;
    
    c=hb-gr * (hb-lb);
    d=lb+gr * (hb-lb);
    
    while abs(c-d)>tol
        
        fc=-func(c,extraIn);
        fd=-func(d,extraIn);
        
        if fc<fd
            hb=d;
            
        else
            lb=c;
        end
        
        c=hb-gr * (hb-lb);
        d=lb+gr * (hb-lb);
    end
    
    np=(hb+lb)/2;
    
end

function [np]=GoldenSection_FindB(lb,hb,tol,M,targ)
    
    gr=(sqrt(5)-1)/2;
    
    c=hb-gr * (hb-lb);
    d=lb+gr * (hb-lb);
    
    while abs(c-d)>tol
        
        fc=abs(targ-CalcTanDel(M,c));
        fd=abs(targ-CalcTanDel(M,d));
        
        if fc<fd
            hb=d;
            
        else
            lb=c;
        end
        
        c=hb-gr * (hb-lb);
        d=lb+gr * (hb-lb);
    end
    
    np=(hb+lb)/2;
    
end

function [loop]=ConstantArea_Klunker(xMin,xMax,A,M,nPoints)
    ratio=xMax-xMin;
    
    Acalc=A;
    
    m = sqrt(M^2-1); %
    Pb = -0.3; % base cp
    
    xl=(1-(m*Pb)/(12*Acalc))/(1-(m*Pb)/(4*Acalc));
    
    t=3/2*Acalc*xl*(1-m*Pb/(12*Acalc));
    
   xCoord=linspace(0,1,ceil(nPoints/2))';
   yCoord=t/2*xCoord/xl.*(2-xCoord/xl);
    
   points=[[xCoord,-yCoord];[xCoord(end:-1:2),yCoord(end:-1:2)]];
   
   ratio=xMax-xMin;
   points(:,1)=points(:,1)*ratio+xMin;
   points(:,2)=points(:,2)/ratio;
   
   [A2]=CalculatePolyArea(points);
   loop.subdivision=points;
    loop.isccw=true;
    loop.subdivision(end+1,:)=loop.subdivision(1,:);
end
