function [] = include_Utilities()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.

    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)

end


%% Parameter Structures Operations

function [trimmedarr]=TrimZeros(arr)
    trimmedarr=arr(arr~=0);
end

function [structdat]=ExploreStructureTree(rootstruct)
    % Recursively explores a structure to return the structures and
    % variables contained within

    maxPosNum=0;
    currentPosVec=[];
    [fields,vars,~]=ExploreStructureTreeRecurse(rootstruct,currentPosVec,maxPosNum);
    structdat.fields=fields;
    structdat.vars=vars;

    function [fields,vars,maxPosNum]=...
            ExploreStructureTreeRecurse...
            (rootstruct,currentPosVec,maxPosNum)
        % Explores a structure tree and returns
        vars=struct([]);
        fields=fieldnames(rootstruct);
        sizStruct=length(fields);
        posNumWorking=maxPosNum;
        maxPosNum=maxPosNum+sizStruct;


        for ii=1:sizStruct
            if isstruct(rootstruct.(fields{ii}))
                workingPosVec=[currentPosVec,posNumWorking+ii]; % Creates the posi
                [newfields,newvars,maxPosNum]=...
                    ExploreStructureTreeRecurse...
                    (rootstruct.(fields{ii})(1),workingPosVec,maxPosNum);
                fields=[fields;newfields];
            else
                workingPosVec=[currentPosVec,posNumWorking+ii];
                newvars.name=fields{ii};
                newvars.vec=workingPosVec;
            end
            vars=[vars,newvars];
            clear newvars
        end


    end

end

function [varargout]=ExtractVariables(varNames,param)

    varInStruct=param.structdat.vardat.names;
    varPosInStruct=param.structdat.vardat.varmatch;
    varargout=cell([1 max(length(varNames),1)]);

    for ii=1:length(varNames)
        %         if strcmp(varNames{ii},'cellLevels')
%             warning('cellLevels used')
%         end
        targLocation=regexp(varInStruct,['#',varNames{ii},'#'], 'once');
        if ~isempty(targLocation)
            targLocation=targLocation+1;
        end
        if isempty(targLocation) || varPosInStruct(targLocation)==0
            targLocation=regexp(varInStruct,varNames{ii}, 'once');
            if isempty(targLocation) || varPosInStruct(targLocation)==0
                error([varNames{ii},' is an invalid variable name'])
            end
        end

        activeVarNum=varPosInStruct(targLocation);
        varargout{ii}=ExtractVariableValue(param,activeVarNum);
    end

    function [varVal]=ExtractVariableValue(param,varNum)

        varstruct=param.structdat.vars(varNum);
        actstruct=param;
        for jj=1:length(varstruct.vec)
            actstruct=actstruct.(param.structdat.fields{varstruct.vec(jj)});
        end
        varVal=actstruct;
    end

end

function [param]=SetVariables(varNames,varValues,param)

    varInStruct=param.structdat.vardat.names;
    varPosInStruct=param.structdat.vardat.varmatch;
    varargout{length(varNames)}=[];

    for ii=1:length(varNames)

        targLocation=regexp(varInStruct,['#',varNames{ii},'#'], 'once');
        if ~isempty(targLocation)
            targLocation=targLocation+1;
        end
        if isempty(targLocation) || varPosInStruct(targLocation)==0
            targLocation=regexp(varInStruct,varNames{ii}, 'once');
            if isempty(targLocation) || varPosInStruct(targLocation)==0
                error([varNames{ii},' is an invalid variable name'])
            end
        end

        activeVarNum=varPosInStruct(targLocation);
        param=SetVariableValue(param,activeVarNum,varValues{ii});
    end

    function [param]=SetVariableValue(param,varNum,newVal)

        varstruct=param.structdat.vars(varNum);
        varPath='param';
        for jj=1:length(varstruct.vec)
            varPath=[varPath,'.',param.structdat.fields{varstruct.vec(jj)}];
        end

        eval([varPath,'=newVal;'])
    end

end

function structdat=GetStructureData(paroptim)

    [structdat]=ExploreStructureTree(paroptim);
    structdat.vardat.names='#';

    for ii=1:length(structdat.vars)
        jj=length(structdat.vardat.names)+1;
        structdat.vardat.names=[structdat.vardat.names,structdat.vars(ii).name,'#'];
        structdat.vardat.varmatch(jj)=ii;
    end

end

%% Vector Operations

function [points]=RemoveIdenticalVectors(points)

    indOrd=1:size(points,1);
    for ii=1:size(points,2)
        [points,iRows]=SortVecColumn(points,ii);
        indOrd=indOrd(iRows);
    end
    [points,indRmv]=RemoveIdenticalConsecutivePoints(points);
    indOrd(indRmv)=[];
    [~,indOrd]=sort(indOrd);
    points=points(indOrd,:);
end

function [vec,iRows]=SortVecColumn(vec,iCol)
    % Sorts according to a columns
    [~,iRows]=sort(vec(:,iCol));
    vec=vec(iRows,:);

end

function [vectorEntries]=RemoveIdenticalEntries(vectorEntries)
    % Function which removes identical entries in a column vector
    vectorEntriesUnsort=vectorEntries;
    [vectorEntries,vectorIndex]=sort(vectorEntries);
    kk=1;
    rmvDI=[];
    for ii=2:length(vectorEntries)
        if vectorEntries(ii)==vectorEntries(ii-1)
            rmvDI(kk)=ii;
            kk=kk+1;
        end
    end
    %vectorEntries(rmvDI)=[];
    vectorIndex(rmvDI)=[];
    vectorEntries=vectorEntriesUnsort(vectorIndex);
end

function [trimmedPoints,indRmv]=RemoveIdenticalConsecutivePoints(points)

%     [~,edgeLength]=LengthProfile(points);
%     indRmv=find(edgeLength<1e-10);
%     indRmv(1)=[]; % remove first point which has 0 distance.
%     trimmedPoints=points;
%     trimmedPoints(indRmv,:)=[];
    indRmv=find((sum((points-points([end,1:end-1],:)).^2,2))<(1e-10)^2);
    if numel(indRmv)==size(points,1);
        indRmv(indRmv==1)=[];
    end
    trimmedPoints=points;
    trimmedPoints(indRmv,:)=[];
end

function cellSimilar=FindIdenticalVector(blockSegments)
    % this function takes in a group of Segments and returns the indices of
    % those identical grouped within a cell array. blockSegments should be a
    % vertical array of horizontal vectors to be compared.

    [m,n]=size(blockSegments);
    blockSegments=sort(blockSegments,2);
    preSortIndex=1:m; % save the index before shuffling
    % Shuffles edges such that similar edges are side by side

    for ii=1:n
        [blockSegments,sortIndex]=SortVecColumn(blockSegments,ii);
        preSortIndex=preSortIndex(sortIndex);
    end
    %compares neighbouring segments
    blockSegTrunc1=blockSegments(1:end-1,:);
    blockSegTrunc2=blockSegments(2:end,:);
    isPrecedent=[0;(sum(blockSegTrunc1==blockSegTrunc2,2)==n)];
    % creates a cell array wih as many elements as there are different edges
    cellSimilar{-sum((isPrecedent-1))}=[];
    kk=0;
    for ii=1:m
        if ~isPrecedent(ii)
            kk=kk+1;
            jj=0;
        end
        jj=jj+1;
        % assigns the presorted index to the similarity array
        cellSimilar{kk}(jj)=preSortIndex(ii);
    end

end

function cellSimilar=FindIdenticalVectorOrd(blockSegments)
    % this function takes in a group of Segments and returns the indices of
    % those identical grouped within a cell array. blockSegments should be a
    % vertical array of horizontal vectors to be compared.

    [m,n]=size(blockSegments);
    preSortIndex=1:m; % save the index before shuffling
    % Shuffles edges such that similar edges are side by side

    for ii=1:n
        [blockSegments,sortIndex]=SortVecColumn(blockSegments,ii);
        preSortIndex=preSortIndex(sortIndex);
    end
    %compares neighbouring segments
    blockSegTrunc1=blockSegments(1:end-1,:);
    blockSegTrunc2=blockSegments(2:end,:);
    isPrecedent=[0;(sum(blockSegTrunc1==blockSegTrunc2,2)==n)];
    % creates a cell array wih as many elements as there are different edges
    cellSimilar{-sum((isPrecedent-1))}=[];
    kk=0;
    for ii=1:m
        if ~isPrecedent(ii)
            kk=kk+1;
            jj=0;
        end
        jj=jj+1;
        % assigns the presorted index to the similarity array
        cellSimilar{kk}(jj)=preSortIndex(ii);
    end

end

%% Various

function [lengthParam,edgeLength]=LengthProfile(points)

    points=points([1,1:end],:);
    pointsVec=points(1:end-1,:)-points(2:end,:);
    edgeLength=sqrt(sum(pointsVec.^2,2));
    lengthParam=cumsum(edgeLength);


end

function [curvParam,edgeCurvNorm,edgeCurve,edgeCurveErr]=CurvatureProfile(points)

    curvFunc=@(pi,pip1,pim1,s1,s2)(-pi.*(s1+s2)+pip1.*s2+pim1.*s1)./(s1.^2.*s2+s2.^2.*s1);
    curvFuncPrec=@(pi,pip1,pim1,s1,s2,eps)(-vpa(pi,eps).*(vpa(s1,eps)+...
        vpa(s2,eps))+vpa(pip1,eps).*vpa(s2,eps)+vpa(pim1,eps).*vpa(s1,eps))./...
        (vpa(s1,eps).^2.*vpa(s2,eps)+vpa(s2,eps).^2.*vpa(s1,eps));

    curvFuncErr=@(pi,pip1,pim1,s1,s2)(-pi.*(s1+s2)+pip1.*s2+pim1.*s1)./(3*s1.*s2+3*s2.*s1);
    normRep=@(v) repmat(sqrt(sum(v.^2,2)),[1 2]);

    edgeCurve=curvFunc(points,points([2:end,1],:),...
        points([end,1:end-1],:),normRep(points-points([2:end,1],:)),...
        normRep(points-points([end,1:end-1],:)));
%     edgeCurve=curvFuncPrec(points,points([2:end,1],:),...
%         points([end,1:end-1],:),normRep(points-points([2:end,1],:)),...
%         normRep(points-points([end,1:end-1],:)),64);
    edgeCurveErr=curvFuncErr(points,points([2:end,1],:),...
        points([end,1:end-1],:),normRep(points-points([2:end,1],:)),...
        normRep(points-points([end,1:end-1],:)));
    edgeCurvNorm=sqrt(sum(edgeCurve.^2,2));
    curvParam=(cumsum(edgeCurvNorm));

end

function [curvParam,edgeCurvNorm,edgeCurve]=CurvatureRadiusProfile(points)

    curvFunc=@(pi,pip1,pim1,s1,s2)(-pi.*(s1+s2)+pip1.*s2+pim1.*s1)./(s1.^2.*s2+s2.^2.*s1);
    normRep=@(v) repmat(sqrt(sum(v.^2,2)),[1 2]);

    edgeCurve=curvFunc(points,points([2:end,1],:),...
        points([end,1:end-1],:),normRep(points-points([2:end,1],:)),...
        normRep(points-points([end,1:end-1],:)));
    edgeCurvNorm=sqrt(sum(edgeCurve.^2,2));
    curvParam=cumsum(edgeCurvNorm);

end


function CopyFileLinux(p1,p2)

    system(['cp -rp ''',p1,''' ''',p2,'''']);
end

function [domainBounds]=MakeCartesianGridBounds(cellLevels)

    %cellLevels=cellLevels+2;
    axRatio=(cellLevels+2)'/cellLevels(1);
    domainBounds=[-axRatio,axRatio];

end

function [domainBounds]=MakeBoundsOuterLayer(cellLevels,domainBounds,isact)

    domainBounds(1,:)= (domainBounds(1,:)-mean(domainBounds(1,:)))*...
        (cellLevels(1)+2+2*(1-isact))/(cellLevels(1)-2*(1-isact))+mean(domainBounds(1,:));
    domainBounds(2,:)= (domainBounds(2,:)-mean(domainBounds(2,:)))*...
        (cellLevels(2)+1)/cellLevels(2)+mean(domainBounds(2,:));
end

%{
function [domainBounds]=MakeCartesianGridBoundsInactE(cellLevels)

    cellNorm=cellLevels;
    cellNorm(1)=cellNorm(1)-2;

    cellLength=1/cellNorm(1);

    axRatio=(cellLevels'+2)*cellLength;
    domainBounds=[-axRatio,axRatio];

end
%}
function [domainBounds]=MakeCartesianGridBoundsInactE(cellLevels)

    cellNorm=cellLevels;
    cellNorm(1)=cellNorm(1)-2;

    cellLength=1/cellNorm(1)/2;

    axRatio=(cellLevels'+2)*cellLength;
    domainBounds=[-axRatio,axRatio];
    domainBounds(1,:)=domainBounds(1,:)+1/2;
end

function [domainBounds]=MakeCartesianGridBoundsInactTE(cellLevels)

    cellNorm=cellLevels;
    cellNorm(1)=cellNorm(1)-1;

    cellLength=1/cellNorm(1)/2;

    axRatio=(cellLevels'+2)*cellLength;
    domainBounds=[-axRatio,axRatio];
    domainBounds(1,:)=domainBounds(1,:)+cellLength+1/2;
end

function [xMin,xMax,t,L,A]=ClosedLoopProperties(points)

    [A]=abs(CalculatePolyArea(points));
    vec=points([end,1:end-1],:)-points;
    L=sum(sqrt(sum(vec.^2,2)));
    t=max(points(:,2))-min(points(:,2));
    xMin=min(points(:,1));
    xMax=max(points(:,1));

end

function [oldIndsNewOrd]=ReverseStructInfo(struct1,fieldShort,fieldLong)
    % this is designed to distribute the content of fieldShort in the same
    % length as the content of fieldLong
   oldIndsNewOrd=cell2mat(cellfun(@(new,old)old*ones([1,numel(new)]),...
                {struct1(:).(fieldLong)},...
                {struct1(:).(fieldShort)},'UniformOutput',false));
end

function [mma]=MovingAverage(dat,span)

    mma=zeros(size(dat));
    mmaspanned=zeros(size(dat));
    for ii=0:span-1
        inds=(1:size(dat,2))-ii;
        indsDest=find(inds>0);
        inds=inds(inds>0);
        mmaspanned(:,indsDest)=mmaspanned(:,indsDest)+1;
        mma(:,indsDest)=mma(:,indsDest)+dat(:,inds);
    end
    mma=mma./mmaspanned;
end

function [mma]=MovingAverageCentral(dat,span)

    mma=zeros(size(dat));
    mmaspanned=zeros(size(dat));
    for ii=-ceil(span/2)+1:ceil(span/2)-1
        inds=(1:size(dat,2))-ii;
        indsDest=find(inds>0 & inds<=size(dat,2));
        inds=inds(inds>0 & inds<=size(dat,2));
        mmaspanned(:,indsDest)=mmaspanned(:,indsDest)+1;
        mma(:,indsDest)=mma(:,indsDest)+dat(:,inds);
    end
    mma=mma./mmaspanned;
end

function [mma]=MovingAverageLoop(dat,span)

    mma=zeros(size(dat));
    mmaspanned=zeros(size(dat));
    for ii=-ceil(span/2)+1:ceil(span/2)-1
        inds=mod((1:size(dat,2))-ii-1,size(dat,2))+1;
        indsDest=find(inds>0 & inds<=size(dat,2));
        inds=inds(inds>0 & inds<=size(dat,2));
        mmaspanned(:,indsDest)=mmaspanned(:,indsDest)+1;
        mma(:,indsDest)=mma(:,indsDest)+dat(:,inds);
    end
    mma=mma./mmaspanned;
end

function [datmat]=MovingIntegralWindowLoop(x,dat,span)

    % calculate start and end of each window
    sizDat=size(dat);
    x=reshape(x,[numel(x),1]);
    dat=reshape(dat,[numel(x),sizDat(sizDat~=numel(x))]);
    xAll=sort(mod([x-span;x;x+span],max(x)));
    datAll = interp1(x,dat,xAll);
    %calculate central portion of integration window
    % the final and first point are assumed to be in the same place
    xallMat=repmat(xAll,[1 numel(x)]);
    xDelta=(xallMat-repmat(x',[numel(xAll) 1]));
    actMat=(mod(xDelta,max(x))>=-span & mod(xDelta,max(x))<span) | ...
        (mod(-xDelta,max(x))>=-span & mod(-xDelta,max(x))<span);
    datmat=repmat(reshape(datAll,[size(datAll,1) 1 size(datAll,2)]),[1 numel(x)]);
    datmat=(datmat+datmat([2:end,1],:,:))/2.*...
        repmat(mod(xallMat([2:end,1],:)-xallMat,max(x)),[1 1 size(dat,2)]);
    datmat(~actMat)=0;
    datmat=reshape(sum(datmat,1),[size(datmat,2) size(datmat,3)])/(2*span);

end

function [datmat]=MovingIntegralWindowLoop2(x,dat,span)

    % calculate start and end of each window
    sizDat=size(dat);
    x=reshape(x,[numel(x),1]);
    dat=reshape(dat,[numel(x),sizDat(sizDat~=numel(x))]);
    xAll=sort(mod([x-span;x;x+span],max(x)));
    datAll = interp1(x,dat,xAll);
    %calculate central portion of integration window
    % the final and first point are assumed to be in the same place
    xallMat=repmat(xAll,[1 numel(x)]);
    xDelta=(xallMat-repmat(x',[numel(xAll) 1]));
    actMat=(mod(xDelta,max(x))>=-span & mod(xDelta,max(x))<span) | ...
        (mod(-xDelta,max(x))>=-span & mod(-xDelta,max(x))<span);
    datmat=repmat(reshape(datAll,[size(datAll,1) 1 size(datAll,2)]),[1 numel(x)]);
%     datmat=(datmat+datmat([2:end,1],:,:))/2.*...
%         repmat(mod(xallMat([2:end,1],:)-xallMat,max(x)),[1 1 size(dat,2)]);
    datmat=(datmat+datmat([2:end,1],:,:))/2;
    datmat(~actMat)=0;
    datmat=reshape(sum(datmat,1),[size(datmat,2) size(datmat,3)])/(2*span);

end

%% From FileExchange




% function [A]=CalculatePolyArea(points)
%
%     pointsVec=points';
%     pointsVec=pointsVec(:);
%     plot(points(:,1),points(:,2));
%     n=length(points(:,1));
%     centreMat=eye(2*n);
%     centreMat=(centreMat+centreMat(:,[end-1:end,1:end-2]))*0.5;
%
%     [rotDif]=[0 -1 0 1; 1 0 -1 0];
%     normMat=zeros(2*n);
%     for ii=1:n-1
%         normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
%     end
%     ii=n;
%     normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
%     normMat((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);
%     A=0.5*(normMat*pointsVec)'*(centreMat*pointsVec);
%
% end

%{
function [A]=CalculatePolyArea(points)

    pointsVec=points';
    pointsVec=pointsVec(:);
    plot(points(:,1),points(:,2));
    n=length(points(:,1));
    centreMat=eye(2*n);
    centreMat=(centreMat+centreMat(:,[end-1:end,1:end-2]))*0.5;

    [rotDif]=[0 -1 0 1; 1 0 -1 0];
    normMat=zeros(2*n);
    for ii=1:n-1
        normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+4))=rotDif;
    end
    ii=n;
    normMat((2*(ii-1)+1):(2*(ii-1)+2),(2*(ii-1)+1):(2*(ii-1)+2))=rotDif(:,1:2);
    normMat((2*(ii-1)+1):(2*(ii-1)+2),1:2)=rotDif(:,3:4);
    A=0.5*(normMat*pointsVec)'*(centreMat*pointsVec);

end
%}
