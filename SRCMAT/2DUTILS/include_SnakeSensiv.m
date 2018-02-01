function [] = include_SnakeSensiv()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end



%% From SnakesSensitivity.m
function [snaxOrd]=CloseRepeatingSnaxLoops(snaxOrd)
    ii=1;
    jj=1;
    while (ii<=numel(snaxOrd))
        
        activeLoop=snaxOrd{ii};
        subActCell=sort(FindObjNum([],activeLoop(jj),activeLoop));
        if numel(subActCell)>1
            snaxOrd{ii}=activeLoop(subActCell(1):subActCell(2)-1);
            snaxOrd{numel(snaxOrd)+1}=activeLoop([subActCell(2):end,1:subActCell(1)-1]);
            jj=0;
        end
        
        ii=ii+floor(jj/length(activeLoop));
        jj=mod(jj,length(activeLoop))+1;
    end
end

function [cellordstruct]=BuildCellConnectivity(snaxel,baseGrid,gridRefined,volfraconnec)
    
    edgeInd=[gridRefined.edge(:).index];
    %activeCellInd=[baseGrid.cell(logical([baseGrid.cell(:).isactive])).index];
    activeCellInd=[baseGrid.cell(:).index];
    [oldCellInd,newCellInd]=NewToOldCell(volfraconnec);
    
    for ii=1:numel(snaxel)
        snaxcellinfo(ii).index=snaxel(ii).index;
        snaxcellinfo(ii).precsnax=snaxel(ii).snaxprec;
        snaxcellinfo(ii).nextsnax=snaxel(ii).snaxnext;
        snaxcellinfo(ii).newcell=gridRefined.edge(...
            FindObjNum([],snaxel(ii).edge,edgeInd)).cellindex;
        
        snaxcellinfo(ii).oldcell=(oldCellInd(...
            FindObjNum([],snaxcellinfo(ii).newcell,newCellInd)));
    end
    
    [snaxOrd]=SplitSnaxLoops(snaxel);
    
    for ii=1:length(snaxOrd)
        cellsRaw=vertcat(snaxcellinfo(snaxOrd{ii}).oldcell);
        cellOrd{ii}=OrderRawCells(cellsRaw,activeCellInd);
    end
    delCell=false(size(cellOrd));
    for ii=1:length(cellOrd)
        delCell(ii)=isempty(cellOrd{ii});
    end
    cellOrd(delCell)=[];
    [cellOrd]=CloseRepeatingSnaxLoops(cellOrd);
    
    cellordstruct=repmat(struct('index',[],'nextcell',[],'prevcell',[],...
        'loopindex',[],'looplength',[]),[1,numel([cellOrd{:}])]);
    kk=1;
    
    for ii=1:numel(cellOrd)
        nCell=numel(cellOrd{ii});
        for jj=1:nCell
            cellordstruct(kk).index=cellOrd{ii}(jj);
            cellordstruct(kk).nextcell=cellOrd{ii}(mod(jj,nCell)+1);
            cellordstruct(kk).prevcell=cellOrd{ii}(mod(jj-2,nCell)+1);
            cellordstruct(kk).loopindex=ii;
            cellordstruct(kk).looplength=nCell;
            
            kk=kk+1;
        end
    end
    
    
end

function [cellsOrdRaw]=OrderRawCells(cellsRaw,activeCellInd)
    
    
    cellsOrdRaw=zeros([1,numel(cellsRaw)]);
    [l2Ind,l1Ind]=find((ones([2,1])*cellsRaw(1,:))==(cellsRaw(2,:)'*ones([1,2]))...
        & (ones([2,1])*cellsRaw(1,:))~=(ones([2,1])*cellsRaw(end,:)));
    if isempty(l2Ind)
        [l2Ind,l1Ind]=find((ones([2,1])*cellsRaw(1,:))==(cellsRaw(2,:)'*ones([1,2])));
    end
    nRows=size(cellsRaw,1);
    cellsOrdRaw(1)=cellsRaw(1,l1Ind(1));
    cellsOrdRaw(2)=cellsRaw(2,l2Ind(1));
    kk=3;
    for ii=2:nRows
        l1Ind=abs(l2Ind(1)-3);
        l2Ind=find(cellsRaw(ii,l1Ind)==cellsRaw(mod(ii,nRows)+1,:));
        
        if isempty(l2Ind)
            error('Non matching cell lists')
        end
        
        cellsOrdRaw(kk)=cellsRaw(ii,l1Ind);
        cellsOrdRaw(kk+1)=cellsRaw(mod(ii,nRows)+1,l2Ind(1));
        kk=kk+2;
        
    end
    cellsOrdRaw=RemoveIdenticalConsecutivePoints(cellsOrdRaw')';
    cellsOrdRaw=cellsOrdRaw(FindObjNum([],cellsOrdRaw,activeCellInd)~=0);
    if cellsOrdRaw(end)==cellsOrdRaw(1)
        cellsOrdRaw(end)=[];
    end
end

function [oldCellInd,newCellInd]=NewToOldCell(volfraconnec)
    
    newCellInd=[volfraconnec.cell(:).newCellInd];
    oldCellInd=ones(size(newCellInd));
    kk=1;
    for ii=1:length(volfraconnec.cell)
        nNew=numel(volfraconnec.cell(ii).newCellInd);
        oldCellInd(kk:kk+nNew-1)=ones([1,nNew])*volfraconnec.cell(ii).oldCellInd;
        kk=kk+nNew;
    end
end

function [snaxOrd]=SplitSnaxLoops(snaxel)
    % Splits a snake into its component loops
    
    kk=1;
    jj=1;
    snaxInd=[snaxel(:).index];
    %snaxOrd=zeros(size(snaxel));
    ordList=1:numel(snaxel);
    ll=1;
    for ii=1:length(snaxel)
        kk=FindObjNum([],[snaxel(kk).snaxnext],snaxInd);
        if ordList(kk)==0
            jj=jj+1;
            kk=min(ordList(ordList~=0));
            ll=1;
        end
        snaxOrd{jj}(ll)=kk;
        ordList(kk)=0;
        ll=ll+1;
    end
    
end
