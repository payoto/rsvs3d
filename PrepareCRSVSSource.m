function []=PrepareCRSVSSource(fileName)
    % THis function generates the C subroutines necessary for the
    % calculation and differentiation of Volume and area in 3-Dimensions
    %
    % Volume contributions are calculated for a triangle built from an
    % polygon
    % Areas are similarly calculated for triangular portions of a polygon
    
    [ccodeGen(5)]=Code_SurfCentroid();
    [ccodeGen(6)]=Code_SurfCentroidDiv();
    [ccodeGen(7)]=Code_SurfCentroidFull();
    [ccodeGen(4)]=Code_SurfIntersect();
    [ccodeGen(1)]=Code_Volume();
    [ccodeGen(2)]=Code_Area();
    [ccodeGen(3)]=Code_LengthEdge();
    for ii=1:numel(ccodeGen)
        ccodeGen(ii)=ChangeUnderScoresToBrackets(ccodeGen(ii));
    end
    
    cellStr=cell(0);
    for ii=1:numel(ccodeGen)
        [cellStr]=[cellStr,CodeGenToCpp(ccodeGen(ii))];
    end
    
    FID=fopen(fileName,'w');
    WriteToFile(cellStr,FID)
    fclose(FID);
    
end
%% Code Generation functions

function ccodeGen=ChangeUnderScoresToBrackets(ccodeGen)
    pRBrack={'([0-9])__','$1]'};
    pLBrack={'_([0-9])','[$1'};
    
    ccodeGen.f=regexprep(ccodeGen.f,pRBrack{1},pRBrack{2});
    ccodeGen.f=regexprep(ccodeGen.f,pLBrack{1},pLBrack{2});
    
    ccodeGen.df=regexprep(ccodeGen.df,pRBrack{1},pRBrack{2});
    ccodeGen.df=regexprep(ccodeGen.df,pLBrack{1},pLBrack{2});
    
    ccodeGen.ddf=regexprep(ccodeGen.ddf,pRBrack{1},pRBrack{2});
    ccodeGen.ddf=regexprep(ccodeGen.ddf,pLBrack{1},pLBrack{2});
end

function [ccodeGen]=CCodeGenStruct(varargin)
    
    [cinputGen]=CInputGenStruct(varargin{:});
    ccodeGen=struct('name','','type','void','return','','f','','df','','ddf','','inputs',cinputGen);
    
end

function [cinputGen]=CInputGenStruct(type,varargin)
    % Generates a structure with the 
    if nargin==0
        type='';
    end
    cinputGen=struct('type',type,'name','');
    cinputGen=repmat(cinputGen,[1 numel(varargin)]);
    
    for ii=1:numel(varargin)
        cinputGen(ii).name=varargin{ii};
    end
    
end

function [cellStr]=CodeGenToCpp(ccodeGen)
    
    fieldsFunc={'f','df','ddf'};
    kk=0;
    for ii=1:numel(fieldsFunc)
        if ~isempty(ccodeGen.(fieldsFunc{ii}))
            kk=kk+1;cellStr{kk}=[ccodeGen.type,' ',ccodeGen.name,'_',fieldsFunc{ii},'('];
            outVar=regexp(ccodeGen.(fieldsFunc{ii}),';\s*\w+','match');
            outVar=deblank(outVar{end}(4:end));
            for jj=1:numel(ccodeGen.inputs)
                cellStr{kk}=[cellStr{kk},ccodeGen.inputs(jj).type,' ',...
                    regexprep(ccodeGen.inputs(jj).name,'#MATCHOUT#',outVar)];
                if jj<numel(ccodeGen.inputs)
                    cellStr{kk}=[cellStr{kk},' , '];
                end
            end
            cellStr{kk}=[cellStr{kk},' ) {'];
            kk=kk+1;cellStr{kk}=ccodeGen.(fieldsFunc{ii});


            kk=kk+1;cellStr{kk}='}';
            kk=kk+1;cellStr{kk}='';
            kk=kk+1;cellStr{kk}='';
        end
    end
    
    
end

function cchar=custcccode(ctarg)
    
    fileName='tempccode.cpp.tmp';
    
    ccode(ctarg,'file',fileName);
    
    FID=fopen(fileName,'r');
    cchar='';
    tempVars=cell(0);
    
    while(~feof(FID))
        chartemp=fgets(FID);
        tempTest=regexp(deblank(chartemp),'^\s*t[1-9][0-9]*','match');
        if ~isempty(tempTest)
            tempVars{end+1}=[tempTest{1},', '];
        end
        
        cchar=[cchar,chartemp];
    end
    if ~isempty(tempVars)
        tempVars=[tempVars{:}];
        cchar=[sprintf('double %s;\n', tempVars(1:end-2)),cchar];
    end
    
    fclose(FID);
    
end

%% Matrix operations for the Geometric properties
% These calculations accept 2 or 3 coordinate triplets and the derivatives
% are provided with regard to these triplets.

% Volume Calculations

function [ccodeGen]=Code_Volume()

    
    p0=sym('p0_%d__',[3,1]);assume(p0,'real');
    p1=sym('p1_%d__',[3,1]);assume(p1,'real');
    p2=sym('p2_%d__',[3,1]);assume(p2,'real');
    
    ccodeGen=CCodeGenStruct('vector<double>','p0','p1','p2');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='Volume';
    
    vol=symfun(p0'*(cross((p2-p0),(p1-p0)))/6,[p0;p1;p2]);
    volJac=jacobian(vol);
    volHes=hessian(vol);
    
    ccodeGen.f=custcccode(vol);
    ccodeGen.df=custcccode(volJac);
    ccodeGen.ddf=custcccode(volHes);
    
end

% Area Calculations


function [ccodeGen]=Code_Area()

    
    p0=sym('p0_%d__',[3,1]);assume(p0,'real');
    p1=sym('p1_%d__',[3,1]);assume(p1,'real');
    p2=sym('p2_%d__',[3,1]);assume(p2,'real');
    
    ccodeGen=CCodeGenStruct('vector<double>','p0','p1','p2');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='Area';
    
    area=symfun(sqrt(transpose(cross((p2-p0),(p1-p0)))*cross((p2-p0),(p1-p0)))/2,[p0;p1;p2]);
    areaJac=jacobian(area);
    areaHes=hessian(area);
    
    ccodeGen.f=custcccode(area);
    ccodeGen.df=custcccode(areaJac);
    ccodeGen.ddf=custcccode(areaHes);
    
end

% Profile Length Calculation


function [ccodeGen]=Code_LengthEdge()

    
    p0=sym('p0_%d__',[3,1]);assume(p0,'real');
    p1=sym('p1_%d__',[3,1]);assume(p1,'real');
    
    ccodeGen=CCodeGenStruct('vector<double>','p0','p1');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='LengthEdge';
    
    lengthedge=symfun(sqrt(sum((p1-p0).^2)),[p0;p1]);
    lengthedgeJac=jacobian(lengthedge);
    lengthedgeHes=hessian(lengthedge);
    
    ccodeGen.f=custcccode(lengthedge);
    ccodeGen.df=custcccode(lengthedgeJac);
    ccodeGen.ddf=custcccode(lengthedgeHes);
    
end

%% Operations for Geometric features

% Surf Object centroid (vec<x>,vec<y>,vec<z>)


function [ccodeGen]=Code_SurfCentroid()

    n=4;
    x=sym('x_%d__',[1 n]);assume(x,'real');
    y=sym('y_%d__',[1 n]);assume(y,'real');
    z=sym('z_%d__',[1 n]);assume(z,'real');
    
    ccodeGen=CCodeGenStruct('vector<double>','x','y','z');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='SurfCentroid';
    
    lengthedge=symfun((sqrt(sum(([x;y;z]-[x([end,1:end-1]);y([end,1:end-1]);z([end,1:end-1])]).^2 ...
        ,1))*transpose([x;y;z])),[x,y,z]);
    lengthedgeJac=jacobian(lengthedge);
    %lengthedgeHes=[hessian(lengthedge(1)),hessian(lengthedge(2)),hessian(lengthedge(3))];
    
    ccodeGen.f=custcccode(lengthedge);
    ccodeGen.df=custcccode(lengthedgeJac);
    ccodeGen.ddf='';
    
end
function [ccodeGen]=Code_SurfCentroidDiv()

    n=4;
    x=sym('x_%d__',[1 n]);assume(x,'real');
    y=sym('y_%d__',[1 n]);assume(y,'real');
    z=sym('z_%d__',[1 n]);assume(z,'real');
    
    ccodeGen=CCodeGenStruct('vector<double>','x','y','z');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='SurfCentroidDiv';
    
    lengthedge=symfun((1)/sum(sqrt(sum(([x;y;z]-[x([end,1:end-1]);...
        y([end,1:end-1]);z([end,1:end-1])]).^2 ,1))),[x,y,z]);
    lengthedgeJac=jacobian(lengthedge);
    %lengthedgeHes=[hessian(lengthedge(1)),hessian(lengthedge(2)),hessian(lengthedge(3))];
    
    ccodeGen.f=custcccode(lengthedge);
    ccodeGen.df=custcccode(lengthedgeJac);
    ccodeGen.ddf='';
    
end


function [ccodeGen]=Code_SurfCentroidFull()

    n=5;
    x=sym('x_%d__',[1 n]);assume(x,'real');
    y=sym('y_%d__',[1 n]);assume(y,'real');
    z=sym('z_%d__',[1 n]);assume(z,'real');
    
    ccodeGen=CCodeGenStruct('vector<double>','x','y','z');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='SurfCentroidDiv';
    
    lengthedge=symfun((sqrt(sum(([x;y;z]-[x([end,1:end-1]);y([end,1:end-1]);z([end,1:end-1])]).^2 ...
        ,1))*transpose([x;y;z]))/sum(sqrt(sum(([x;y;z]-[x([end,1:end-1]);...
        y([end,1:end-1]);z([end,1:end-1])]).^2 ,1))),[x([1:2,4]),y([1:2,4]),z([1:2,4])]);
    lengthedgeJac=jacobian(lengthedge);
    %lengthedgeHes=[hessian(lengthedge(1)),hessian(lengthedge(2)),hessian(lengthedge(3))];
    
    ccodeGen.f=custcccode(lengthedge);
    ccodeGen.df=custcccode(lengthedgeJac);
    ccodeGen.ddf='';
    
end

% snaxedge to triangular edge intersection



function [ccodeGen]=Code_SurfIntersect()

    
     
    p1=sym('p1_%d__',[3,1]);assume(p1,'real');
    p2=sym('p2_%d__',[3,1]);assume(p2,'real');
    v0=sym('v0_%d__',[3,1]);assume(v0,'real');
    v01=sym('v01_%d__',[3,1]);assume(v01,'real');
    v02=sym('v02_%d__',[3,1]);assume(v02,'real');
    
    v11=sym('v11_%d__',[3,1]);assume(v11,'real');
    v12=sym('v12_%d__',[3,1]);assume(v12,'real');
    
    c=sym('c_%d__',[3,1]);assume(c,'real');
    
    ccodeGen=CCodeGenStruct('vector<double>','p1','p2','v0','v01','v02','v11','v12','c');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='SurfIntersect';
    
    dg1=(v11-v01);
    dg2=(v12-v02);
    
    nv=cross(dg1,dg2);%opt1
    nv=cross((p2-p1),dg1+dg2);%opt2
    
    ns=cross(nv,(p2-p1));
    
    lengthedge=symfun(((dot(ns,p1))-(dot(ns,c)))/(dot(ns,v0-c)),[p1;p2]);
    lengthedgeJac=jacobian(lengthedge);
    lengthedgeHes=hessian(lengthedge);
    
    ccodeGen.f=custcccode(lengthedge);
    ccodeGen.df=custcccode(lengthedgeJac);
    ccodeGen.ddf=custcccode(lengthedgeHes);
    
end

