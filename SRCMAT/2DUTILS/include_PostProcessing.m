function [] = include_PostProcessing()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end

%% General
function []=savefig(figh,figname)

print(figh,figname,'-djpeg','-r600')

end

function [marker,t]=GenerateResultMarker(typDat)
    
    t=now;
    marker=[datestr(t,'yyyy-mm-ddTHHMMSS')...
        ,'_',typDat];
    
end

function [resultDirectory]=GenerateResultDirectoryName(marker,resultRoot,...
        archiveName,t)
    if ~exist('t','var'),t=now;end
    dateSubFolders=['Archive_',datestr(now,'yyyy_mm'),'\Day_',datestr(t,29)];
    resultDirectory=[resultRoot,filesep,archiveName,filesep,dateSubFolders,...
        filesep,'Dir_',marker];
    
    resultDirectory=MakePathCompliant(resultDirectory);
    
    mkdir(resultDirectory)
end

function pathName=MakePathCompliant(pathName)
    
    compStr=computer;
    if strcmp(compStr(1:2),'PC')
        pathName=regexprep(pathName,'/','\\');
    else
        
        pathName=regexprep(pathName,'\\','/');
    end
end

function []=WriteToFile(cellLoops,FID)
    % writes the data to the file
    
    
    for ii=1:length(cellLoops)
        if numel(cellLoops{ii})>0
            for jj=1:length(cellLoops{ii}(:,1))
                fprintf(FID,'%s \n',cellLoops{ii}(jj,:));
            end
        else
            fprintf(FID,'\n');
        end
    end
end

function []=ReorganiseSubPlots(ax,sizSubPlot,outerPad,interPad,fontSizes)
    
    if ~exist('outerPad','var'),outerPad=[0.1,0.1,0.05,0.05];end
    if ~exist('interPad','var'),interPad=[0.05,0.05];end
    if ~exist('fontSizes','var'),fontSizes=[10,12];end
    
    xHeight=1-(outerPad(1)+outerPad(2));
    xLength=(xHeight-interPad(1)*(sizSubPlot(1)-1))/sizSubPlot(1);
    xStart=outerPad(1);
    
    yHeight=1-(outerPad(3)+outerPad(4));
    yLength=(yHeight-interPad(2)*(sizSubPlot(2)-1))/sizSubPlot(2);
    yStart=1-outerPad(4)-yLength;
    
    for ii=1:length(ax)
        jj=ii;
        xInd=rem(jj-1,sizSubPlot(1));
        yInd=ceil(jj/sizSubPlot(1))-1;
        
        loXCorn=(xInd)*(xLength+interPad(1))+xStart;
        loYCorn=-(yInd)*(yLength+interPad(2))+yStart;
        
        posVec=[loXCorn,loYCorn,xLength,yLength];
        set(ax(ii),'position',posVec)
        set(ax(ii),'fontsize',fontSizes(1))
        set(get(ax(ii),'xlabel'),'fontsize',fontSizes(2))
        set(get(ax(ii),'ylabel'),'fontsize',fontSizes(2))
        set(get(ax(ii),'zlabel'),'fontsize',fontSizes(2))
        
    end
    
end

function []=CreateValidFolder(pathName)
    error('replace this shit by mkdir')
    
%     if strcmp(compStr(1:2),'PC')
%         pathName=regexprep(pathName,'/','\\');
%         system(['md "',pathName,'"']);
%     else
%         
%         pathName=regexprep(pathName,'\\','/');
%         system(['mkdir ''',pathName,''''])
%     end
end

%% Plotting 

function [datCol]=ProjectColormap(cMap,cDat,cBounds)
    
    if nargin==2
        cBounds=[min(cDat),max(cDat)];
    end
    if cBounds(1)==cBounds(2)
        cBounds(2)=1+cBounds(2);
    end
    nCol=length(cMap(:,1));
    
    cDat(cDat>max(cBounds))=max(cBounds);
    cDat(cDat<min(cBounds))=min(cBounds);
    
    datMap=linspace(cBounds(1),cBounds(2),nCol)';
    
    
    datCol=interp1(datMap,cMap,cDat);
    
end

%% Parameter File

function []=GenerateParameterFile(FID,param,t,marker)
    
    paramCell{1}='% Parameter File';
    paramCell{2}=['% ',datestr(t)];
    paramCell{3}=['% ',marker];
    for ii=1:length(param.structdat.vars)
        paramCell=[paramCell,ExtractVariablePathAndValue2(param,ii)];
        
    end
    
    WriteToFile(paramCell,FID);
    fclose(FID);
end

function [paramStr]=ExtractVariablePathAndValue(param,varNum)
    
    varstruct=param.structdat.vars(varNum);
    actstruct=param;
    pathVar='param';
    for jj=1:length(varstruct.vec)
        pathVar=[pathVar,'.',param.structdat.fields{varstruct.vec(jj)}];
        actstruct=actstruct.(param.structdat.fields{varstruct.vec(jj)});
    end
    
    varVal=actstruct;
    [varStr]=GenerateVariableString(varVal);
    paramStr=[pathVar,' = ',varStr,';'];
end

function [paramStr]=ExtractVariablePathAndValue2(param,varNum)
    
    varstruct=param.structdat.vars(varNum);
    actstruct=param;
    pathVar{1}='param';
    
    
    for jj=1:length(varstruct.vec)
        tempPath=[];
        tempPath{numel(pathVar)}=[];
        for kk=1:numel(pathVar)
            pathVar{kk}=[pathVar{kk},'.',param.structdat.fields{varstruct.vec(jj)}];
            nTestVar=numel(eval(pathVar{kk}));
            if (jj<length(varstruct.vec)) && nTestVar>1
                for ii=1:nTestVar
                    tempPath{kk}{ii}=[pathVar{kk},'(',int2str(ii),')'];
                end
            else
                tempPath{kk}{1}=pathVar{kk};
            end
        end
        pathVar=[tempPath{:}];
    end
    paramStr{numel(pathVar)}=[];
    for ii=1:numel(pathVar)
        [varStr]=GenerateVariableString(eval(pathVar{ii}));
        paramStr{ii}=[pathVar{ii},' = ',varStr,';'];
    end
end

function [varStr]=GenerateVariableString(startVar)
    
    classVar=class(startVar);
    [m,n]=size(startVar);
    varStr='';
    switch classVar
        case 'char'
            
            openStr='[';
            closeStr=']';
            if ~isempty(startVar)
                for ii=1:m
                    varStrCell{ii,1}=['''',startVar(ii,:),''''];
                end
                [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,1);
                else
                varStr='''''';
            end
        case 'cell'
            
            openStr='{';
            closeStr='}';
            for ii=1:m
                for jj=1:n
                    varStrCell{ii,jj}=GenerateVariableString(startVar{ii,jj});
                end
            end
            if exist('varStrCell','var')
                [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
            else
                [varStr]=[openStr,closeStr];
            end
        case 'double'
            
            openStr='[';
            closeStr=']';
            [varStrCell{1:max([m 1]),1:max([n 1])}]=deal(' ');
            for ii=1:m
                for jj=1:n
                    varStrCell{ii,jj}=num2str(startVar(ii,jj),24);
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
        case 'logical'
            openStr='[';
            closeStr=']';
            for ii=1:m
                for jj=1:n
                    if startVar(ii,jj) 
                        curStr='true';
                    else
                        curStr='false';
                    end
                    varStrCell{ii,jj}=curStr;
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
        otherwise
            if ~isempty(regexp(classVar,'int','once'))
                
            openStr='[';
            closeStr=']';
            for ii=1:m
                for jj=1:n
                    varStrCell{ii,jj}=int2str(startVar(ii,jj));
                end
            end
            [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n);
            end
            warning('Class is not catered for and will not be printed correctly to parameter file')
    end
    
end

function [varStr]=RecursiveStringGeneration(openStr,closeStr,varStrCell,m,n)

    
    varStr=openStr;
    if m==0
        m=1;
    end
    if n==0
        n=1;
    end
    for ii=1:m-1
        for jj=1:n-1
            varStr=[varStr,varStrCell{ii,jj},','];
        end
        varStr=[varStr,varStrCell{ii,n},';'];
    end
    for jj=1:n-1
        varStr=[varStr,varStrCell{m,jj},','];
    end
    
    varStr=[varStr,varStrCell{m,n},closeStr];
end

function []=CopyDiary(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    fileName=['diary_',marker,'.log'];
    originalLayFile=[cd,'\Result_Template\Latest_Diary.log'];
    originalLayFile=MakePathCompliant(originalLayFile);
    copyfile(originalLayFile,MakePathCompliant([writeDirectory,filesep,fileName]))
    
end



%% Boundary Output to .dat file

function []=BoundaryOutput(loop,FID,typeLoop,buildInternal)
    
    if nargin<4
        buildInternal=false;
        if nargin<3
            typeLoop='subdivision';
        end
    end
    if numel(loop)==0
        error('Invalid loop input, check ''fill''')
    end
    % trim loops and extract data
    if ~buildInternal
        [isInternal]=FindInternalLoop(loop,typeLoop);
        loop=loop(~isInternal);
    end
    if numel(loop)==0
        error('Loop is Empty after internal Trimming, check ''fill''')
    end
    loopout=TrimLoops(loop,typeLoop);
    % format numeric data to printable strings
    cellLoops=DataToString(loopout);
    
    % print string data to file
    WriteToFile(cellLoops,FID)
    
end

function [G]=BuildMatlabGeometryMatrix(loop,typeLoop)
    [isInternal]=FindInternalLoop(loop,typeLoop);
    
    nPts=0;
    for ii=1:numel(loop)
        nPts=nPts+size(loop(ii).(typeLoop),1);
    end
    
    G=zeros(7,nPts);
    nDone=0;
    for ii=1:numel(loop)
        nDoneNew=nDone+size(loop(ii).(typeLoop),1);
        padOnes=ones([1,size(loop(ii).(typeLoop),1)]);
        actCoord=loop(ii).(typeLoop);
        if ~loop(ii).isccw
            actCoord=flip(actCoord);
        end
        actCoord=actCoord';
        G(:,nDone+1:nDoneNew)=[padOnes*2;
            actCoord(1,:);
            actCoord(1,[2:end,1]);
            actCoord(2,:);
            actCoord(2,[2:end,1]);
            padOnes*mod(isInternal(ii)+1,2);
            padOnes*mod(isInternal(ii),2)];
        nDone=nDoneNew;
    end
    
end

function []=OutputFreefempp(loop,typeLoop)
    [isInternal]=FindInternalLoop(loop,typeLoop);
    
    
    
end

function [loopout]=TrimLoops(loop,typeLoop)
    % function extracting the data that must be written to the boundary.dat
    % file
    if nargin<2
        typeLoop='subdivision';
    end
    
    nLoop=length(loop);
    
    for ii=1:nLoop
        endInd=length(loop(ii).(typeLoop)(:,1));
        
        isNeg2=sum(sum(loop(ii).(typeLoop)(1:2,:)==loop(ii).(typeLoop)(end-1:end,:)))...
            ==numel(loop(ii).(typeLoop)(1:2,:));
        isNeg1=sum(loop(ii).(typeLoop)(1,:)==loop(ii).(typeLoop)(end,:))...
            ==numel(loop(ii).(typeLoop)(1,:));
        
        if isNeg2
            endInd=endInd-2;
        elseif isNeg1
            endInd=endInd-1;
        else
        
        end
        
        if loop(ii).isccw
            loopout.surf(ii).coord=loop(ii).(typeLoop)(1:endInd,:);
        else
            loopout.surf(ii).coord=loop(ii).(typeLoop)(endInd:-1:1,:);
        end
        loopout.surf(ii).nvertex=size(loopout.surf(ii).coord,1);
        loopout.surf(ii).nfaces=loopout.surf(ii).nvertex;
    end
    
    loopout.total.nvertex=sum([loopout.surf(:).nvertex]);
    loopout.total.nfaces=sum([loopout.surf(:).nfaces]);
    loopout.total.nloops=nLoop;
    
end

function cellLoops=DataToString(loopout)
    % Transforms the data organised in loopout into a cell array of strings
    % ready to be written to a file
    
    cellLoops{1}=int2str(loopout.total.nloops);
    cellLoops{2}=int2str([loopout.total.nvertex, loopout.total.nfaces]);
    kk=3;
    
    % run through the different surfaces
    for ii=1:loopout.total.nloops
        cellLoops{kk}=int2str(loopout.surf(ii).nvertex);
        kk=kk+1;
        % Write an array of numbers in the correct format
        for jj=1:length(loopout.surf(ii).coord(:,1))
            cellLoops{kk}='';
            for ll=1:length(loopout.surf(ii).coord(jj,:))
                cellLoops{kk}=[cellLoops{kk},...
                    num2str(loopout.surf(ii).coord(jj,ll),24),'  '];
            end
            kk=kk+1;
        end
        
    end
    
end

function [isInternal]=FindInternalLoop(loop,typeLoop)
    % this function works by checkking all the inequalities around a
    % polygon
    isInternal=zeros(size(loop));
    for ii=1:length(loop)
        polygonPoints=loop(ii).(typeLoop);
%         polyVec=polygonPoints([2:end,1],:)-polygonPoints;
%         polyNorm=[polyVec(:,2),-polyVec(:,1)];
%         
%         numCond=size(polygonPoints,1);
        for jj=[1:ii-1,ii+1:length(loop)]
            
%             comparePoint=ones([size(polyNorm,1),1])*loop(jj).subdivision(1,:);
%             conditionNum=sum(polyNorm.*(comparePoint-polygonPoints),2);
%             allPos=sum(conditionNum>=0)==numCond;
%             allNeg=sum(conditionNum<=0)==numCond;
              isIn=inpolygon(loop(jj).(typeLoop)(1,1),loop(jj).(typeLoop)(1,2),...
                  loop(ii).(typeLoop)(:,1),loop(ii).(typeLoop)(:,2));
            isInternal(jj)=isInternal(jj) + (isIn);
        end
    end
    
end

function [loop]=BoundaryInput(fileName)
    
    fid=fopen(fileName,'r');
    
    if ~feof(fid)
        numLoops=str2double(fgetl(fid));
        numPoints=str2num(fgetl(fid));
        
        loop=repmat(struct('coord',[],'isccw',true,'isinternal',false),[1,numLoops]);
        
        for ii=1:numLoops
            numPts(ii)=str2double(fgetl(fid));
            loop(ii).coord=zeros(numPts(ii),2);
            for jj=1:numPts(ii)
                loop(ii).coord(jj,:)=str2num(fgetl(fid));
            end
            
        end
        
        if sum(numPts)~=numPoints(1)
            error('Reading loop in did not pass point sum check')
        end
        
        if ~feof(fid)
            error('Reading loop in did not reach end of file')
        end
    else
        error(['File failed to open: ',fileName])
    end
    fclose(fid);
    
end

%% Displacements Output

function []=DisplacementsOutput(catloop,FIDsurf,FIDdisp,typeLoop,buildInternal)
    
    if nargin<5
        buildInternal=false;
        if nargin<4
            typeLoop='subdivision';
        end
    end
    if numel(catloop)==0
        error('Invalid loop input, check ''fill''')
    end
    % trim loops and extract data
    if ~buildInternal
        [isInternal]=FindInternalLoop(catloop(:,1),typeLoop);
        catloop=catloop(find(~isInternal),:);
    end
    if numel(catloop)==0
        error('Loop is Empty after internal Trimming, check ''fill''')
    end
    [trimpoints,trimdisp]=TrimDisplacementPoints(catloop,typeLoop);
    % format numeric data to printable strings
    cellPoints=DispToString(trimpoints);
    cellDisplacement=DispToString(trimdisp);
    
    % print string data to file
    WriteToFile(cellPoints,FIDsurf)
    WriteToFile(cellDisplacement,FIDdisp)
    
end

function [trimpoints,trimdisp]=TrimDisplacementPoints(catloop,typeLoop)
    % Removes double points and sums them 
    points=vertcat(catloop(:,1).(typeLoop));
    displacements=vertcat(catloop(:,3).(typeLoop));
    
    cellSimilar=FindIdenticalVector(points);
    trimdisp=zeros([numel(cellSimilar),size(displacements,2)]);
    trimpoints=zeros([numel(cellSimilar),size(displacements,2)]);
    for ii=1:numel(cellSimilar)
        trimdisp(ii,:)=sum(displacements(cellSimilar{ii},:),1);
        trimpoints(ii,:)=points(cellSimilar{ii},:);
    end
    
    if size(trimdisp,2)<3
        trimdisp(:,3)=0;
        trimpoints(:,3)=0;
    end
end

function cellLoops=DispToString(array)
    % Transforms the data organised in loopout into a cell array of strings
    % ready to be written to a file
    
    
    kk=1;
    cellLoops{kk}=int2str(size(array,1));
    kk=kk+1;
    % Write an array of numbers in the correct format
    for jj=1:size(array,1)
        cellLoops{kk}='';
        for ll=1:length(array(jj,:))
            cellLoops{kk}=[cellLoops{kk},...
                num2str(array(jj,ll),24),'  '];
        end
        kk=kk+1;
    end
    
    
    
end
%% Video Output functions

function []=MakeVideo(movStruct,fps,quality,fileName)
    
    writerObj = VideoWriter(fileName);
    writerObj.FrameRate=fps;
    writerObj.Quality=quality;
    open(writerObj)
    writeVideo(writerObj,movStruct)
    close(writerObj)
end

%% Layout Operation functions

function []=PersnaliseLayFile(FID,pltFile)
    
    frewind(FID);
    layData{1}='#!MC 1410';
    layData{2}=['$!VarSet |LFDSFN1| = ''"',pltFile,'"'''];
    if iscell(pltFile)
        for ii=1:length(pltFile)
            layData{1+ii}=['$!VarSet |LFDSFN',int2str(ii),'| = ''"',pltFile{ii},'"'''];
        end
    end
    WriteToFile(layData,FID)
end

%% File Opening Functions

function [FID]=OpenBoundaryFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['boundary_',marker,'.dat'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
end

function [FID,fileName]=OpenTecPLTFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['tec360dat_',marker,'.plt'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
end

function [FID]=OpenTecLayFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['tec360lay_',marker,'.lay'];
    originalLayFile=[cd,'\Result_Template\Layout_Template.lay'];
    originalLayFile=MakePathCompliant(originalLayFile);
    copyfile(originalLayFile,[writeDirectory,filesep,fileName])
    FID=fopen([writeDirectory,filesep,fileName],'r+');
    
end

function [FID]=OpenParamFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['param_',marker,'.dat'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
end

function [FID]=OpenCommentsFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['Comments_',marker,'.txt'];
    FID=fopen([writeDirectory,filesep,fileName],'w+');
    
end

function [FID]=OpenIndexFile(resultRoot,archiveName)
    % Creates a file in the current directory to write data to.
    
    resultRoot=MakePathCompliant(resultRoot);
    fileName=['Index_',archiveName,'.txt'];
    FID=fopen([resultRoot,filesep,archiveName,filesep,fileName],'a');
    
end

function [FID]=NameVideoFile(writeDirectory,marker)
    % Creates a file in the current directory to write data to.
    
    writeDirectory=MakePathCompliant(writeDirectory);
    fileName=['Video_',marker,'.avi'];
    FID=[writeDirectory,filesep,fileName];
    
end

%% Comments File

function [indexEntry]=MakeCommentsFile(FID,param,t,resultDirectory)
    
    varExtract={'typDat','case','noteFiles','tags'};
    [typDat,caseStr,noteFiles,tags]=ExtractVariables(varExtract,param);
    
    
    headerLines=GenerateCommentHeader(t,resultDirectory,typDat,caseStr,tags);
    automatedComments=ConcatenateAutomaticComments(noteFiles);
    indexEntry{1}=GenerateIndexEntry(t,resultDirectory,typDat,caseStr,tags);
    
    kk=1;
    breakLines{kk}=['--------------------------------------------------'];
    kk=kk+1;
    breakLines{kk}=['ADDITIONAL COMMENTS:'];
    kk=kk+1;
    breakLines{kk}=[' '];
    
    WriteToFile(headerLines,FID)
    WriteToFile(automatedComments,FID)
    WriteToFile(breakLines,FID)
    
end

function automatedComments=ConcatenateAutomaticComments(noteFiles)
    
    for ii=length(noteFiles):-1:1
        noteFileName{ii}=MakePathCompliant([cd,'\Result_Template\Notes_',noteFiles{ii},'.txt']);
        fidNote=fopen(noteFileName{ii},'r');
        rawComments(ii)=textscan(fidNote,'%s','Delimiter','\n');
        fclose(fidNote);
    end
    kk=1;
    automatedComments{kk}=['--------------------------------------------------'];
    kk=kk+1;
    automatedComments{kk}=['PRESET COMMENTS:'];
    kk=kk+1;
    automatedComments{kk}=[' '];
    for ii=1:length(noteFiles)
        kk=kk+1;
        automatedComments{kk}=['-------'];
        kk=kk+1;
        automatedComments{kk}=noteFileName{ii};
        kk=kk+1;
        automatedComments{kk}=['-------'];
        kk=kk+1;
        kk2=kk+length(rawComments{ii})-1;
        automatedComments(kk:kk2)=rawComments{ii};
        kk=kk2;
        kk=kk+1;
        automatedComments{kk}=[' '];
        
    end
    kk=kk+1;
    automatedComments{kk}=[' '];
    
end

function headerLines=GenerateCommentHeader(t,resultDirectory,typDat,caseStr,tags)
    
    kk=1;
    headerLines{kk}=datestr(t);
    kk=kk+1;
    headerLines{kk}=resultDirectory;
    kk=kk+1;
    headerLines{kk}=['Data File: ',typDat];
    kk=kk+1;
    headerLines{kk}=['Case: ',caseStr];
    kk=kk+1;
    headerLines{kk}=[' '];
    kk=kk+1;
    headerLines{kk}=['Tags: '];
    for ii=1:length(tags)
        headerLines{kk}=[headerLines{kk},tags{ii},', '];
    end
    kk=kk+1;
    headerLines{kk}=[' '];
    
    
end

function indexLine=GenerateIndexEntry(t,resultDirectory,typDat,caseStr,tags)
    
    indexLine='';
    indexLine=[indexLine,datestr(t)];
    indexLine=[indexLine,', '];
    indexLine=[indexLine,typDat];
    indexLine=[indexLine,', '];
    indexLine=[indexLine,caseStr];
    indexLine=[indexLine,', '];
    for ii=1:length(tags)-1
        indexLine=[indexLine,tags{ii},' - '];
    end
    indexLine=[indexLine,tags{ii+1}];
    indexLine=[indexLine,', '];
    indexLine=[indexLine,resultDirectory];
    
end

function []=WriteFullInfoProfile(writeDirectory,nProf,marker,t,nIter)
    fid=fopen([writeDirectory,filesep,'InfoProfile_',int2str(nIter),'_',int2str(nProf),'.dat'],'w');
    
    ii=1;
    cellDat{ii}=marker;ii=ii+1;
    cellDat{ii}=datestr(t);ii=ii+1;
    cellDat{ii}=['iteration : ',int2str(nIter)]; ii=ii+1;
    cellDat{ii}=['Profile : ',int2str(nProf)]; ii=ii+1;
    
    WriteToFile(cellDat,fid);
    fclose(fid);
end

%% Error Report File

function [FID]=OpenErrorReportFile(rootOptim,marker)
    % Creates a file in the current directory to write data to.
    
    rootOptim=MakePathCompliant(rootOptim);
    fileName=['ErrorReport_',marker,'.txt'];
    FID=fopen([rootOptim,filesep,fileName],'a');
    
end

function []=GenerateErrorReportFile(t,marker,FID)
    
    paramCell{1}='# Error Report File';
    paramCell{2}=['# ',datestr(t)];
    paramCell{3}=['# ',marker];
    paramCell{4}=[' '];
    
    WriteToFile(paramCell,FID);
    fclose(FID);
end

function []=GenerateErrorReportEntries(fID,nIter,errorReports,indexEntries)
    
    writeReport{1}='------------------------------------------------------------';
    writeReport{2}=['   ITERATION ',int2str(nIter)];
    writeReport{3}='------------------------------------------------------------';
    kk=4;
    for ii=1:length(errorReports)
        
        if ~isempty(errorReports{ii})
            writeReport{kk}=indexEntries{ii};
            kk=kk+1;
            writeReport{kk}=errorReports{ii};
            kk=kk+1;
        end
        
    end

    WriteToFile(writeReport,fID);
    fclose(fID);
    
end

%% Generate Restart Binary

function []=GenerateRestartBinary(resultDirectory,marker,restartstruct)
    
    fileName=[resultDirectory,filesep,'restart_',marker,'.mat'];
    save(fileName,'-struct','restartstruct');
    
end

function []=GenerateProfileBinary(resultDirectory,marker,restartstruct)
    
    fileName=[resultDirectory,filesep,'restart_',marker,'.mat'];
    fileName=MakePathCompliant(fileName);
    save(fileName,'-struct','restartstruct');
    
end

%% Tecplot

function [cellMesh]=CellEdgeMesh(coordDat,vertIndex,connDat,vectorDat,strandID,time)
    
    
    [numNode,nDim]=size(coordDat);
    [numElm,nConn]=size(connDat);
    % Zone Header
    if ~exist('time','var')
        cellHeader=FELINESEGHeader(numNode,numElm,strandID);
    else
        cellHeader=FELINESEGHeader(numNode,numElm,strandID,time);
    end
    % Coordinates
    for ii=numNode:-1:1
        cellNodeCoord{ii}='';
        for jj=1:nDim
            cellNodeCoord{ii}=[cellNodeCoord{ii},...
                num2str(coordDat(ii,jj),24),'  '];
        end
        for jj=1:nDim+1
            cellNodeCoord{ii}=[cellNodeCoord{ii},...
                num2str(vectorDat(ii,jj),24),'  '];
        end
    end
    % Connectivity
    for ii=numElm:-1:1
        cellConnCoord{ii}='';
        for jj=1:nConn
            cellConnCoord{ii}=[cellConnCoord{ii},...
                int2str(FindObjNum([],connDat(ii,jj),vertIndex)),'  '];
        end
        
    end
    
    cellMesh=[cellHeader,cellNodeCoord,cellConnCoord];
    
end

%% Tecplot Headers for various zone types

function cellHeader=FELINESEGHeader(numNodes,numElm,strandID,time)
    
%     if ~exist('numNodes','var');numNodes=1;end
%     if ~exist('numElm','var');numElm=1;end
    
    
    nodesStr=['NODES=',int2str(numNodes),' '];
    elmStr=['ELEMENTS=',int2str(numElm)];
    idStr=['STRANDID=',int2str(strandID)];
    
    kk=1;
    
    cellHeader{kk}=['VARIABLES = "X" ,"Y" ,"U" ,"V", "MAG"']; kk=kk+1;
    cellHeader{kk}=['ZONE']; kk=kk+1;
    cellHeader{kk}=['VARLOCATION=([1-5]=NODAL)']; kk=kk+1;
    cellHeader{kk}=nodesStr; kk=kk+1;
    cellHeader{kk}=elmStr; kk=kk+1;
    if exist('time','var');
        cellHeader{kk}=idStr; kk=kk+1;
        timeStr=['SOLUTIONTIME=',num2str(time,'%.24f'),' '];
        cellHeader{kk}=timeStr; kk=kk+1;
    end
    cellHeader{kk}='DATAPACKING=POINT'; kk=kk+1;
    cellHeader{kk}='ZONETYPE=FELINESEG'; kk=kk+1;
    
    
    
end

function cellHeader=CellCentredDataHeader(numNodes,numElm,nFaces,vararg)
    % vararg is a structure that can be used to specify arguments:
    % Possible vararg fields
    % strandID, time, connecShareZone, varsharezone
    
    varFields=fieldnames(vararg);
    for ii=1:length(varFields)
        eval([varFields{ii},'=vararg.(varFields{ii});'])
    end
    
    
    nodesStr=['NODES=',int2str(numNodes),' '];
    elmStr=['ELEMENTS=',int2str(numElm)];
    faceStr=['FACES=',int2str(nFaces)];

    kk=1;
    
    cellHeader{kk}=['VARIABLES = "X" ,"Y" ,"TARGFILL" ,"VOLFRAC", "DIFF"']; kk=kk+1;
    cellHeader{kk}=['ZONE']; kk=kk+1;
    cellHeader{kk}=['VARLOCATION=([1-2]=NODAL ,[3-5]=CELLCENTERED)']; kk=kk+1;
    cellHeader{kk}=nodesStr; kk=kk+1;
    cellHeader{kk}=elmStr; kk=kk+1;
    cellHeader{kk}=faceStr; kk=kk+1;
    cellHeader{kk}=['NUMCONNECTEDBOUNDARYFACES=0']; kk=kk+1;
    cellHeader{kk}=['TOTALNUMBOUNDARYCONNECTIONS=0']; kk=kk+1;
    if ~isempty(time);
        idStr=['STRANDID=',int2str(strandID)];
        cellHeader{kk}=idStr; kk=kk+1;
        timeStr=['SOLUTIONTIME=',num2str(time,'%.24f'),' '];
        cellHeader{kk}=timeStr; kk=kk+1;
    end
    if ~isempty(connecShareZone);
        connecZoneStr=['CONNECTIVITYSHAREZONE=',int2str(connecShareZone),' '];
        cellHeader{kk}=connecZoneStr; kk=kk+1;
    end
    if ~isempty(varsharezone);
        cellHeader{kk}=VariableShareZoneString(varsharezone); kk=kk+1;
    end
    cellHeader{kk}='DATAPACKING=BLOCK'; kk=kk+1;
    cellHeader{kk}='ZONETYPE=FEPOLYGON'; kk=kk+1;
    
end

function [str]=VariableShareZoneString(varsharezone)
    % Creates teh string for Variable sharing depending on the structure
    % varsharezone with fields .vars and .zone
    
    str='VARSHARELIST=(';
    for ii=1:length(varsharezone)
        str=[str,'['];
        for jj=1:length(varsharezone(ii).vars)
            str=[str,int2str(varsharezone(ii).vars(jj))];
            if jj~=length(varsharezone(ii).vars)
                str=[str,','];
            end
        end
        str=[str,']=',int2str(varsharezone(ii).zone)];
        if ii~=length(varsharezone)
            str=[str,','];
        end
    end
    str=[str,')'];
end

%% Incomplete processing

function [outinfo]=ReconstructOutinfo(optimstruct)
    
    allRootDir=repmat(struct('rootDir',''),[1 numel(optimstruct)]);
    rmR=[];
    for ii=1:numel(optimstruct)
        try
            kk=1;
            while isempty(optimstruct(ii).population(kk).location)
                kk=kk+1;
            end
            allRootDir(ii).rootDir=regexprep(optimstruct(ii).population(kk).location,'iteration.*$','');
        catch
            
        
            rmR=[ii];
        end
    end
    allRootDir(rmR)=[];
    [allRootDir,uniqueRootDir]=IdentifyUniqueOptions(allRootDir);
    outinfo=repmat(struct('marker','','tOutput',[],'rootDir',''),[1 numel(uniqueRootDir{1})]);
    [outinfo(:).rootDir]=deal(uniqueRootDir{1}{:});
    for ii=1:numel(outinfo)
        delInd(ii)=isempty(outinfo(ii).rootDir);
    end
    outinfo(delInd)=[];
    for ii=1:numel(outinfo)
        [outinfo(ii).rootDir,outinfo(ii).tOutput,outinfo(ii).marker]=FindTime2(outinfo(ii).rootDir);
    end
    [~,sortOrd]=sort([outinfo(:).tOutput]);
    outinfo=outinfo(sortOrd);
end

function [refinestruct,patternList]=IdentifyUniqueOptions(refinestruct,jobfields)
    
    if nargin<2
        jobfields=fieldnames(refinestruct);
    end
    
    %
    
    n=numel(refinestruct);
    
    patternList=cell(0);
    for ii=1:numel(jobfields);
        patternList{ii}={};
        
        for jj=1:n
            optfound=false;
            
            for kk=1:numel(patternList{ii})
                test=false;
                if ischar(refinestruct(jj).(jobfields{ii}))
                    test=strcmp(patternList{ii}{kk},refinestruct(jj).(jobfields{ii}));
                elseif isnumeric(refinestruct(jj).(jobfields{ii})) || islogical(refinestruct(jj).(jobfields{ii}))
                    test=(patternList{ii}{kk}==refinestruct(jj).(jobfields{ii}));
                end
                if test
                    refinestruct(jj).([jobfields{ii},'num'])=kk;
                    optfound=true;
                    break
                end
            end
            if ~optfound
                patternList{ii}{end+1}=refinestruct(jj).(jobfields{ii});
                refinestruct(jj).([jobfields{ii},'num'])=numel(patternList{ii});
            end
        end
    end
    
    
    
    
end

function [t,marker]=FindTime(pathStr)
    
    [returnPath,returnName]=FindDir(pathStr,'Comments',false);
    
    fid=fopen(returnPath{1},'r');
    
    t=datenum(fgetl(fid));
    
    returnName=regexprep(returnName{1},'Comments_','');
    marker=regexprep(returnName,'.txt','');
    
end

function [pathStr,t,marker]=FindTime2(pathStr)
    
    splitPath=regexp(pathStr,'[\\,/]','split');
    pathStr=regexprep(pathStr,'[\\,/]$','');
    
    dirLoc=flip(find(~cellfun(@isempty,regexp(splitPath,'Dir'))));
    dirName=splitPath{dirLoc};
    splitDir=regexp(dirName,'_');
    t=datenum(dirName(splitDir(1)+1:splitDir(2)-1),'yyyy-mm-ddTHHMMSS');
    marker=dirName(splitDir(2)+1:end);
end

function [paramoptim]=ReconstructParameter(pathStr,marker)
    try
        [returnPath,returnName]=FindDir(pathStr,'FinalParam',0);
        load(returnPath{1})
        if ~exist('paramoptim','var')
            error('Failure to load')
        end
    catch
        
        [rootParam,nameParam]=FindDir(pathStr,'param',0);
        isOptParam=cellfun(@isempty,regexp(nameParam,'_parametrisation'));
        matName=matlab.lang.makeValidName(regexprep(nameParam{isOptParam},'.dat',''));
        copyfile(rootParam{isOptParam},['.',filesep,matName,'.m']);
        eval(matName)
        paramoptim=param;
        paramoptim.structdat=GetStructureData(paramoptim);
        
        copyfile(rootParam{~isOptParam},['.',filesep,matName,'.m']);
        eval(matName)
        paramoptim.parametrisation=param;
        paramoptim.parametrisation.structdat=...
            GetStructureData(paramoptim.parametrisation);
    end
end

function [iterstruct]=ReconstructIterationStructure(pathStr,nIter,paramoptim)
    
    varExtract={'direction','optimMethod'};
    [direction,optimMethod]=ExtractVariables(varExtract,paramoptim);
    
    if numel(nIter)>1
        nIterStart=nIter(1);
        nIterEnd=nIter(2);
        nIter=nIter(2)-nIter(1)+1;
    end
    
    iterstruct=struct([]);
    iterstruct=repmat(iterstruct,[nIter,1]);
    
    for ii=nIterStart:nIterEnd
        
        [iterPath,~]=FindDir(pathStr,['iteration_',int2str(ii)],true);
        
        datName=['population_iteration_',int2str(ii)];
        [returnPath,~]=FindDir(iterPath{1},datName,false);
        load(returnPath{1},'population');
        iterstruct(ii-nIterStart+1).population=population;
        
    end
    
    if TestGreed(optimMethod)
        for ii=2:nIter
            iterPrecedent=[iterstruct(ii-1).population(:).objective];
            iterCurr=[iterstruct(ii).population(:).objective];
            switch direction
                case 'min'
                    iterCopy=iterCurr>iterPrecedent;
                case 'max'
                    iterCopy=iterCurr<iterPrecedent;
            end

            [iterstruct(ii).population(iterCopy)]=deal(iterstruct(ii-1).population(iterCopy));


        end
    end
    
end

function [isgreedy]=TestGreed(optimMethod)
    
   isgreedy=true;
   switch optimMethod
       case 'DE'
       case 'DEtan'
       case 'conjgrad'
           isgreedy=false;
       case 'conjgradls'
           isgreedy=false;
       otherwise
           error('Need to say if optimisation method is greedy or not')
   end
    
    
    
end

%% Triangles support

function [polystruct,structmesh]=OutputLoop2TrianglePoly(fileName,loop,typeLoop,structmesh,nElmZone)
    % Loop is a set of closed loops
    % typeLoop is the field from that structure to use
    % structmesh is a structure with the following fields or a char
    % .boundaryMarkers is used to mark boundaries as internal or external
    %   depending on the final use of the mesh
    %   It is structure =struct('loop',0,'outbound',0,'sym',0) or a char
    % .flagInOut is used to differentiate between meshing inside geometries
    %   and outside (Flow vs structure) 0=structural mesh 1=flow mesh
    % .domainBound is the distance in chords from which the domain should be
    %   generated
    %  Information about .poly file : https://www.cs.cmu.edu/~quake/triangle.poly.html
    if nargin<5
        nElmZone=4000;
    end
    if ischar(structmesh)
        [structmesh]=BoundaryMarkersCharCases(structmesh,nElmZone);
    end
    boundaryMarkers=structmesh.boundaryMarkers;
    flagInOut=structmesh.flagInOut;
    domainBound=structmesh.domainBound;
    
    polystruct=struct('vertex',zeros([0,4]),'segment',zeros([0,4]),'hole',zeros([0,3]),'region',zeros([0,4]));
    % Builds Loops as .poly
    [polystruct]=BuildPolyStruct(polystruct,loop,typeLoop,boundaryMarkers.loop,flagInOut);
    
    % Builds domain boundary if needed
    if flagInOut
        coord=vertcat(loop.(typeLoop));
        dist=max((max(coord)-min(coord)))*domainBound(1);
        ldom.coord=dist*[-1 -1; -1 1; 1 1; 1 -1]+repmat((max(coord)+min(coord))/2,[4,1]);
        [polystruct]=BuildPolyStruct(polystruct,ldom,'coord',...
             boundaryMarkers.outbound,0);
         [p]=FindInternalPoint(ldom.coord);
        tempHole=[size(polystruct.region,1)+1,p,...
            abs(CalculatePolyArea(ldom.coord)/domainBound(2))];
        polystruct.region=[polystruct.region;tempHole];
    end
    
    [polystruct]=TriangleConstructSubDomains(loop,typeLoop,polystruct,structmesh);
    % output to file
    FID=fopen(fileName,'w');
    polyCell=WritePolyStruct2Cell(polystruct);
    WriteToFile(polyCell,FID)
    fclose(FID);
end

function [polystruct]=TriangleConstructSubDomains(loop,typeLoop,polystruct,structmesh)
    
    normRep=@(v) repmat(sqrt(sum(v.^2,2)),[1 2]);
    plotPoints= @(points) points;%plot(points([1:end,1],1),points([1:end,1],2),'-o');
    
    coord=vertcat(loop.(typeLoop));
    ldom.coordbase=coord(convhull(coord(:,1),coord(:,2)),:);
    
    plotPoints(ldom.coordbase);
    ldom.coordbase(1:2:size(ldom.coordbase,1)*2,:)=ldom.coordbase;
    ldom.coordbase(2:2:end+1,:)=(ldom.coordbase(1:2:end,:)+ldom.coordbase([3:2:end,1],:))/2;
     
    plotPoints(ldom.coordbase);
    
    
    ldom.normVec=(ldom.coordbase([2:end,1],:)-ldom.coordbase)./...
        normRep(ldom.coordbase([2:end,1],:)-ldom.coordbase)+...
        (ldom.coordbase-ldom.coordbase([end,1:end-1],:))./...
        normRep(ldom.coordbase-ldom.coordbase([end,1:end-1],:));
    ldom.normVec=([0 1;-1 0]*ldom.normVec')'./normRep(ldom.normVec);
    
    plotPoints(ldom.coordbase);
    dist=[max(ldom.coordbase);min(ldom.coordbase)];
    dist=dist(1,:)-dist(2,:);
    ldom.coord=ldom.coordbase+ldom.normVec*max(dist)*0.05/2;
    plotPoints(ldom.coordbase);
    paramspline.splineCase='convhulltri';
    [ldom.coordbase]=ResampleSpline(ldom.coordbase,paramspline);
    plotPoints(ldom.coordbase);
    ldom.normVec=(ldom.coordbase([2:end,1],:)-ldom.coordbase)./...
        normRep(ldom.coordbase([2:end,1],:)-ldom.coordbase)+...
        (ldom.coordbase-ldom.coordbase([end,1:end-1],:))./...
        normRep(ldom.coordbase-ldom.coordbase([end,1:end-1],:));
    ldom.normVec=([0 1;-1 0]*ldom.normVec')'./normRep(ldom.normVec);
    
    
    for ii=1:size(structmesh.subDomains,1) %find(isfinite(structmesh.subDomains(:,1)))'%
        m=structmesh.subDomains(ii,1);
        if isfinite(m) % build a new region and hole
            ldom.coord=ldom.coordbase+ldom.normVec*max(dist)*m/2;
             plotPoints(ldom.coord);


            if structmesh.subDomains(ii,2)>0
                Acons=abs(CalculatePolyArea(ldom.coord)/structmesh.subDomains(ii,2));
            elseif structmesh.subDomains(ii,2)<0
                Acons=-structmesh.subDomains(ii,2);
            end
            [~,edgeLength]=LengthProfile(ldom.coord);

            paramspline.samplingN=ceil(min([sum(edgeLength)/sqrt(Acons),100])/2)*2+1;

            [ldom.coord]=ResampleSpline(ldom.coord,paramspline);
            plotPoints(ldom.coord);
            [polystruct]=BuildPolyStruct(polystruct,ldom,'coord',0,0);
            [p]=FindInternalPoint(ldom.coord);
        else % if it is infinite put a subdomain constraint on the other side
            [~,p]=FindInternalPoint(ldom.coord);
        end
        if structmesh.subDomains(ii,2)>0
            tempHole=[size(polystruct.region,1)+1,p,...
                abs(CalculatePolyArea(ldom.coord)/structmesh.subDomains(ii,2))];
        elseif structmesh.subDomains(ii,2)<0
            tempHole=[size(polystruct.region,1)+1,p,...
               -structmesh.subDomains(ii,2)];
        end
            
        polystruct.region=[polystruct.region;tempHole];
    end
    
    %         ldom.coord=ldom.coordbase-[-1 -1; -1 1; 1 1; 1 -1]*max(dist)*m/2;
%         ldom.coord(:,1)=(ldom.coord(:,1)-mean(ldom.coord(:,1)))*...
%             (max(r,1)*min(m,2)+max(m-2,2))*max(r,1)...
%             +mean(ldom.coord(:,1));
%         ldom.coord(:,2)=(ldom.coord(:,2)-mean(ldom.coord(:,2)))*...
%             (max(1/r,1)*min(m,2)+max(m-2,2))*max(1/r,1)...
%             +mean(ldom.coord(:,2));
end

function [structmesh2]=BoundaryMarkersCharCases(structmesh,nElmZone)
    switch structmesh
        case 'cutcellflow'
            structmesh2.boundaryMarkers=struct('loop',-1,'outbound',-2,'sym',-3);
            structmesh2.flagInOut=1;
            structmesh2.domainBound=[25 10];
            % each line of subdomains are is a distance cage and the number
            % of cells that must fit.
            %             structmesh2.subDomains=[.1 2*nElmZone;.25 nElmZone;.5 nElmZone;...
            %                 1 nElmZone;2 nElmZone;8 nElmZone];
            structmesh2.subDomains=logspace(log10(.1),log10(20),15)';
            if nElmZone>0
                structmesh2.subDomains(:,2)=logspace(log10(2),log10(1),15)'*nElmZone;
            elseif nElmZone<0
                structmesh2.subDomains(:,2)=logspace(log10(1),...
                    log10(structmesh2.subDomains(end,1)^2),15)'*nElmZone;
            end
            
        case 'cutcellflow2'
            structmesh2.boundaryMarkers=struct('loop',-1,'outbound',-2,'sym',-3);
            structmesh2.flagInOut=1;
            structmesh2.domainBound=[25 10];
            % each line of subdomains are is a distance cage and the number
            % of cells that must fit.
            %             structmesh2.subDomains=[.1 2*nElmZone;.25 nElmZone;.5 nElmZone;...
            %                 1 nElmZone;2 nElmZone;8 nElmZone];
            nRef=round((log2(max(structmesh2.domainBound)^2)-log2(abs(nElmZone)))/2);
            structmesh2.subDomains(1:nRef,2)=nElmZone*4.^(0:nRef-1)';
            
            structmesh2.subDomains(:,1)=2*sqrt(abs(structmesh2.subDomains(1:nRef,2))).*logspace(log10(20),...
                log10(10),nRef)';
            outOffBound=find(structmesh2.subDomains(:,1)>max(structmesh2.domainBound));
            structmesh2.subDomains(outOffBound(1),1)=-inf;
            structmesh2.subDomains(outOffBound(2:end),:)=[];
            
        case 'su2'
            
            error('Not coded yet')
            case 'structures'
                error('Not coded yet')
        case 'none'
            structmesh2.boundaryMarkers=struct('loop',0,'outbound',0,'sym',0);
            structmesh2.flagInOut=1;
            structmesh2.domainBound=25;
            structmesh2.subDomains=zeros(0,2);
        otherwise
            warning('unrecognised boundaryMarkers')
            boundaryMarkers=struct('loop',1,'outbound',1,'sym',1);
    end
    
    
end

function [polyCell]=WritePolyStruct2Cell(polystruct)
    
    ii=1;
    
    dim=2;
    isBoundMarker=1;
    %  <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
    polyCell{ii}=sprintf('%i %i %i %i',size(polystruct.vertex,1),dim,...
        size(polystruct.vertex,2)-(dim+1+isBoundMarker),isBoundMarker);
    ii=ii+1;
    polyCell{ii}=num2str(polystruct.vertex,'%i  %.24e    %.24e    %i');
    ii=ii+1;
    polyCell{ii}=sprintf('%i  %i',size(polystruct.segment,1),isBoundMarker);
    ii=ii+1;
    polyCell{ii}=num2str(polystruct.segment,'%i  %i    %i   %i');
    ii=ii+1;
    polyCell{ii}=sprintf('%i',size(polystruct.hole,1));
    ii=ii+1;
    polyCell{ii}=num2str(polystruct.hole,'%i  %.24e    %.24e');
    ii=ii+1;
    polyCell{ii}=sprintf('%i',size(polystruct.region,1));
    ii=ii+1;
    polyCell{ii}=num2str(polystruct.region,'%i  %.24e    %.24e   %.24e');
    
end

function [polystruct]=BuildPolyStruct(polystruct,loop,typeLoop,loopBoundMark,flagInOut)
    % currently no support for "Attributes"
    % First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
    % Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
    % One line: <# of segments> <# of boundary markers (0 or 1)>
    % Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
    % One line: <# of holes>
    % Following lines: <hole #> <x> <y>
    % Optional line: <# of regional attributes and/or area constraints>
    % Optional following lines: <region #> <x> <y> <attribute> <maximum area>
    [isInternal]=FindInternalLoop(loop,typeLoop);
    isInternal=xor(isInternal,flagInOut);
    
    for ii=1:numel(loop)
        tempVert=loop(ii).(typeLoop);
        
        tempVert=[size(polystruct.vertex,1)+[1:size(tempVert,1)]',...
            tempVert,ones([size(tempVert,1),1])*loopBoundMark]; %#ok<AGROW>
        
        tempSegment=[size(polystruct.segment,1)+[1:size(tempVert,1)]',...
            [tempVert(:,1),tempVert([2:end,1],1)]...
            ,ones([size(tempVert,1),1])*loopBoundMark];
        polystruct.vertex=[polystruct.vertex;tempVert];
        polystruct.segment=[polystruct.segment;tempSegment];
        if isInternal(ii)
            [p]=FindInternalPoint(loop(ii).(typeLoop));
            tempHole=[size(polystruct.hole,1)+1,p];
            polystruct.hole=[polystruct.hole;tempHole];
        end
    end

end

function [pIn,pOut]=FindInternalPoint(coord)
    np=size(coord,1);
    eps=1e-6;
    epsNorm=1e-12;
    flagDone=false;
    while ~flagDone
        for jj=1:np
            p1=coord(jj,:);
            p2=coord(mod(jj,np)+1,:);
            pr=(p1+p2)/2;
            pvec=p1-p2;
            if sqrt(sum(pvec.^2))>epsNorm
                
                pvec=([0 1; -1 0]*(pvec/sqrt(sum(pvec.^2))*eps)')';
                pTest=repmat(pr,[2 1])+[1;-1]*pvec;
                isIn=inpolygon(pTest(:,1),pTest(:,2),coord(:,1),coord(:,2));
                if any(isIn)
                    flagDone=true;
                    break;
                end
            end
        end
        if ~epsNorm
            flagDone=true;
        end
        if ~flagDone
            eps=1e-10;
            epsNorm=0;
        end
    end
    pIn=pTest(find(isIn),:);
    pOut=pTest(find(~isIn),:);
    if isempty(pIn)
        errstruct.message='Could not find an internal point to specify the hole of .poly file';
        errstruct.identifier='mesh:triangle:polyFile:hole';
        error(errstruct)
    end
end