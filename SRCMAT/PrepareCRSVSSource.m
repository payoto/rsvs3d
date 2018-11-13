function []=PrepareCRSVSSource(targetPath,caseName)
    % THis function generates the C subroutines necessary for the
    % calculation and differentiation of Volume and area in 3-Dimensions
    %
    % Volume contributions are calculated for a triangle built from an
    % polygon
    % Areas are similarly calculated for triangular portions of a polygon
    if true
        [ccodeGen(1)]=Code_Volume_d();
        [ccodeGen(2)]=Code_Volume();
%         [ccodeGen(2)]=Code_Area();
%         [ccodeGen(3)]=Code_LengthEdge();
%         [ccodeGen(4)]=Code_SurfIntersect();
%         [ccodeGen(5)]=Code_SurfCentroid4();
%         [ccodeGen(6)]=Code_SurfCentroid5();
%         [ccodeGen(7)]=Code_SurfCentroid6();
%         [ccodeGen(8)]=Code_SurfCentroidConnec();
%         [ccodeGen(9)]=Code_SurfCentroidNoConnec();
%         [ccodeGen(10)]=Code_SurfCentroidSelf();
    else
        [ccodeGen(8)]=Code_SurfCentroidConnec();
        [ccodeGen(9)]=Code_SurfCentroidNoConnec();
        [ccodeGen(10)]=Code_SurfCentroidSelf();
    end
    
    for ii=1:numel(ccodeGen)
        ccodeGen(ii)=ChangeUnderScoresToBrackets(ccodeGen(ii));
    end
    
    [cellSource,cellHeader]=IncludeLib(caseName);
    
    for ii=1:numel(ccodeGen)
        [cellSourceTemp,cellHeaderTemp]=CodeGenToCpp(ccodeGen(ii));
        [cellSource]=[cellSource,cellSourceTemp];
        [cellHeader]=[cellHeader,cellHeaderTemp];
    end
    cellSource=regexprep(cellSource,'double &   A0','ArrayVec<double> &   A0');
    cellHeader=regexprep(cellHeader,'double &   A0','ArrayVec<double> &   A0');
    cellHeader{end+1}='#endif';
    
    sourceName=[targetPath,filesep,caseName,'.cpp'];
    headerName=[targetPath,filesep,caseName,'.hpp'];
    
    FID=fopen(sourceName,'w');
    WriteToFile(cellSource,FID)
    fclose(FID);
    FID=fopen(headerName,'w');
    WriteToFile(cellHeader,FID)
    fclose(FID);
    for ii=1:numel(ccodeGen)
        for jj=fieldnames(ccodeGen(ii).symfun)'
            MakeTexEquation(ccodeGen(ii).symfun.(jj{1}),['Code_',ccodeGen(ii).name,'_',jj{1},'.tex'],...
                'C:\Users\ap1949\Local Documents\PhD\Write Up\Snakes_Documentation\texmath')
        end
    end
    
end

%% Set pre-code
function [cellSource,cellHeader]=IncludeLib(caseName)
    
    cellSource={'#include <vector>','#include <cmath>','#include "vectorarray.hpp"',...
        'using namespace std;',''};
    
    cellHeader=[{['#ifndef ',upper(caseName),'_H_INCLUDED'],...
        ['#define ',upper(caseName),'_H_INCLUDED'],''},cellSource];
    cellSource=[cellSource,{'#include "RSVSmath_automatic.hpp"'}];
    
end


%% Code Generation functions

function ccodeGen=ChangeUnderScoresToBrackets(ccodeGen)
    pRBrack={'([0-9]*)__','$1-1]'};
    pLBrack={'_([0-9]*)','[$1'};
    
    ccodeGen.f=regexprep(ccodeGen.f,pRBrack{1},pRBrack{2});
    ccodeGen.f=regexprep(ccodeGen.f,pLBrack{1},pLBrack{2});
    [ccodeGen.f]=ChangeMinus1toNum(ccodeGen.f);
    
    ccodeGen.df=regexprep(ccodeGen.df,pRBrack{1},pRBrack{2});
    ccodeGen.df=regexprep(ccodeGen.df,pLBrack{1},pLBrack{2});
    [ccodeGen.df]=ChangeMinus1toNum(ccodeGen.df);
    
    ccodeGen.ddf=regexprep(ccodeGen.ddf,pRBrack{1},pRBrack{2});
    ccodeGen.ddf=regexprep(ccodeGen.ddf,pLBrack{1},pLBrack{2});
    [ccodeGen.ddf]=ChangeMinus1toNum(ccodeGen.ddf);
end

function [str]=ChangeMinus1toNum(str)
    
    [a,b,c]=regexp(str,'\[[0-9]-1\]','start','end','match');
    
    for ii=1:numel(a)
       c{ii}=['[',int2str(str2double(c{ii}(2:end-3))-1),']'];
       str(a(ii):b(ii))=[c{ii},blanks(1+b(ii)-a(ii)-numel(c{ii}))];
    end
    
end

function [ccodeGen]=CCodeGenStruct(varargin)
    
    [cinputGen]=CInputGenStruct(varargin{:});
    ccodeGen=struct('name','','type','void','return','','f','','df','','ddf','',...
        'inputs',cinputGen,'symfun',struct([]));
    
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

function [cellSource,cellHeader]=CodeGenToCpp(ccodeGen)
    
    fieldsFunc={'f','df','ddf'};
    kk=0;
    for ii=1:numel(fieldsFunc)
        if ~isempty(ccodeGen.(fieldsFunc{ii}))
            kk=kk+1;cellSource{kk}=[ccodeGen.type,' ',ccodeGen.name,'_',fieldsFunc{ii},'('];
            outVar=regexp(ccodeGen.(fieldsFunc{ii}),';\s*\w+','match');
            outVar=deblank(outVar{end}(4:end));
            for jj=1:numel(ccodeGen.inputs)
                cellSource{kk}=[cellSource{kk},ccodeGen.inputs(jj).type,' ',...
                    regexprep(ccodeGen.inputs(jj).name,'#MATCHOUT#',outVar)];
                if jj<numel(ccodeGen.inputs)
                    cellSource{kk}=[cellSource{kk},' , '];
                end
            end
            cellHeader{ii}=[cellSource{kk},' );'];
            cellSource{kk}=[cellSource{kk},' ) {'];
            kk=kk+1;cellSource{kk}=ccodeGen.(fieldsFunc{ii});


            kk=kk+1;cellSource{kk}='}';
            kk=kk+1;cellSource{kk}='';
            kk=kk+1;cellSource{kk}='';
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
    
    ccodeGen=CCodeGenStruct('const vector<double>&','p0','p1','p2');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='Volume';
    
    vol=symfun(p0'*(cross((p2-p0),(p1-p0)))/6  ...
        ,[p0;p1;p2]);
%         + p1'*(cross((p0-p1),(p2-p1)))/6 + ...
%         p2'*(cross((p1-p2),(p0-p2)))/6 ...
    volJac=jacobian(vol);
    volHes=hessian(vol);
    
    ccodeGen.f=custcccode(vol);
    ccodeGen.df=custcccode(volJac);
    ccodeGen.ddf=custcccode(volHes);
    ccodeGen.symfun(1).fun=vol;
    ccodeGen.symfun.jac=volJac;
    ccodeGen.symfun.hes=volHes;
end

function [ccodeGen]=Code_Volume_d()

    syms d0 d1 d2;
    g0s=sym('g0s_%d__',[3,1]);
    g1s=sym('g1s_%d__',[3,1]);
    g2s=sym('g2s_%d__',[3,1]);
    g0e=sym('g0e_%d__',[3,1]);
    g1e=sym('g1e_%d__',[3,1]);
    g2e=sym('g2e_%d__',[3,1]);
    assume(d0,'real');
    assume(d1,'real');
    assume(d2,'real');
    assume(g0e,'real');
    assume(g1e,'real');
    assume(g2e,'real');
    assume(g0s,'real');
    assume(g1s,'real');
    assume(g2s,'real');
    
    p0=symfun(d0*(g0e-g0s)+g0s,[d0;g0s;g0e;d1;g1s;g1e;d2;g2s;g2e]);
    p1=symfun(d1*(g1e-g1s)+g1s,[d0;g0s;g0e;d1;g1s;g1e;d2;g2s;g2e]);
    p2=symfun(d2*(g2e-g2s)+g2s,[d0;g0s;g0e;d1;g1s;g1e;d2;g2s;g2e]);
    ccodeGen=CCodeGenStruct('const vector<double>&','g0s','g1s','g2s','g0e','g1e','g2e');
    [ccodeGen.inputs]=[CInputGenStruct('double ','d0', 'd1', 'd2'), ...
        ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='Volume2';
    
    vol=symfun(p0'*(cross((p2-p0),(p1-p0)))/6  ...
        ,[d0;d1;d2]);
%         + p1'*(cross((p0-p1),(p2-p1)))/6 + ...
%         p2'*(cross((p1-p2),(p0-p2)))/6 ...
    volJac=jacobian(vol);
    volHes=hessian(vol);
    
    ccodeGen.f=custcccode(vol);
    ccodeGen.df=custcccode(volJac);
    ccodeGen.ddf=custcccode(volHes);
    ccodeGen.symfun(1).fun=vol;
    ccodeGen.symfun.jac=volJac;
    ccodeGen.symfun.hes=volHes;
end

% Area Calculations


function [ccodeGen]=Code_Area()

    
    p0=sym('p0_%d__',[3,1]);assume(p0,'real');
    p1=sym('p1_%d__',[3,1]);assume(p1,'real');
    p2=sym('p2_%d__',[3,1]);assume(p2,'real');
    eps=sym('eps'); assume(eps,'real'); assume(eps>0);
    ccodeGen=CCodeGenStruct('const vector<double>&','p0','p1','p2');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='Area';
    
    area=symfun(sqrt(transpose(cross((p2-p0),(p1-p0)))*cross((p2-p0),(p1-p0))+eps)/2,[p0;p1;p2]);
    areaJac=jacobian(area);
    areaHes=hessian(area);
    
    ccodeGen.f=custcccode(area);
    ccodeGen.df=custcccode(areaJac);
    ccodeGen.ddf=custcccode(areaHes);
    ccodeGen.symfun(1).fun=area;
    ccodeGen.symfun.jac=areaJac;
    ccodeGen.symfun.hes=areaHes;
    
end

% Profile Length Calculation


function [ccodeGen]=Code_LengthEdge()

    
    p0=sym('p0_%d__',[3,1]);assume(p0,'real');
    p1=sym('p1_%d__',[3,1]);assume(p1,'real');
    
    ccodeGen=CCodeGenStruct('const vector<double>&','p0','p1');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='LengthEdge';
    
    lengthedge=symfun(sqrt(sum((p1-p0).^2)),[p0;p1]);
    lengthedgeJac=jacobian(lengthedge);
    lengthedgeHes=hessian(lengthedge);
    
    ccodeGen.f=custcccode(lengthedge);
    ccodeGen.df=custcccode(lengthedgeJac);
    ccodeGen.ddf=custcccode(lengthedgeHes);
    ccodeGen.symfun(1).fun=lengthedge;
    ccodeGen.symfun.jac=lengthedgeJac;
    ccodeGen.symfun.hes=lengthedgeHes;
    
end

%% Operations for Geometric features

% Surf Object centroid (vec<x>,vec<y>,vec<z>)

function [ccodeGen]=Code_SurfCentroidSelf()

    n=3;
    [ccodeGen]=Code_SurfCentroid(n,[2]);
    ccodeGen.name='SurfCentroidSelf';
    
end
function [ccodeGen]=Code_SurfCentroid4()

    n=4;
    [ccodeGen]=Code_SurfCentroid(n,1:n);
    ccodeGen.name=['SurfCentroid',int2str(n)];
    
end
function [ccodeGen]=Code_SurfCentroid5()

    n=5;
    [ccodeGen]=Code_SurfCentroid(n,1:n);
    ccodeGen.name=['SurfCentroid',int2str(n)];
    
end
function [ccodeGen]=Code_SurfCentroid6()

    n=6;
    [ccodeGen]=Code_SurfCentroid(n,1:n);
    ccodeGen.name=['SurfCentroid',int2str(n)];
    
end

function [ccodeGen]=Code_SurfCentroidDiv()

    n=5;
    x=sym('x_%d__',[1 n]);assume(x,'real');
    y=sym('y_%d__',[1 n]);assume(y,'real');
    z=sym('z_%d__',[1 n]);assume(z,'real');
    
    ccodeGen=CCodeGenStruct('const vector<double>&','x','y','z');
    [ccodeGen.inputs]=[ccodeGen.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='SurfCentroidDiv';
    
    lengthedge=symfun((1)/sum(sqrt(sum(([x;y;z]-[x([end,1:end-1]);...
        y([end,1:end-1]);z([end,1:end-1])]).^2 ,1))),[x,y,z]);
    lengthedgeJac=jacobian(lengthedge);
    %lengthedgeHes=[hessian(lengthedge(1)),hessian(lengthedge(2)),hessian(lengthedge(3))];
    
    ccodeGen.f=custcccode(lengthedge);
    ccodeGen.df=custcccode(lengthedgeJac);
    ccodeGen.ddf='';
    
    ccodeGen.symfun(1).fun=lengthedge;
    ccodeGen.symfun.jac=lengthedgeJac;
    %ccodeGen.symfun.hes=;
    
end


function [ccodeGen]=Code_SurfCentroidAll()

    n=5;
    [ccodeGen]=Code_SurfCentroid(n,[1 2 4]);
    ccodeGen.name='SurfCentroidFull';
    
end

function [ccodeGen]=Code_SurfCentroidConnec()

    n=6;
    [ccodeGen]=Code_SurfCentroid(n,[2 3]);
    ccodeGen.name='SurfCentroidConnec';
    
end
function [ccodeGen]=Code_SurfCentroidNoConnec()

    n=8;
    [ccodeGen]=Code_SurfCentroid(n,[2 6]);
    ccodeGen.name='SurfCentroidNoConnec';
    
    
end

function [ccodeGen]=Code_SurfCentroidNoConnec2()

    n=20;
    [ccodeGen]=Code_SurfCentroid(n,[10 11]);
    ccodeGen.name='SurfCentroidNoConnec';
    
    
end

function [ccodeGen]=Code_SurfCentroid(n,ind)

    x=sym('x_%d__',[1 n]);assume(x,'real');
    y=sym('y_%d__',[1 n]);assume(y,'real');
    z=sym('z_%d__',[1 n]);assume(z,'real');
    
    ccodeGen=CCodeGenStruct('const vector<double>&','x','y','z');
    ccodeGen2=CCodeGenStruct('const double','totD','X_dot_d','Y_dot_d','Z_dot_d');
    [ccodeGen.inputs]=[ccodeGen.inputs,ccodeGen2.inputs,CInputGenStruct('double &','#MATCHOUT#')];
    ccodeGen.name='SurfCentroid';
    
    varArray=[x(ind),y(ind),z(ind)];
    
    X_dot_d=(sqrt(sum(([x;y;z]-[x([end,1:end-1]);y([end,1:end-1]);z([end,1:end-1])]).^2 ...
        ,1))*transpose(([x]+[x([end,1:end-1])])/2));
    Y_dot_d=(sqrt(sum(([x;y;z]-[x([end,1:end-1]);y([end,1:end-1]);z([end,1:end-1])]).^2 ...
        ,1))*transpose(([y]+[y([end,1:end-1])])/2));
    Z_dot_d=(sqrt(sum(([x;y;z]-[x([end,1:end-1]);y([end,1:end-1]);z([end,1:end-1])]).^2 ...
        ,1))*transpose(([z]+[z([end,1:end-1])])/2));
    totD=sum(sqrt(sum(([x;y;z]-[x([end,1:end-1]);...
        y([end,1:end-1]);z([end,1:end-1])]).^2 ,1)));
    
    crep=@(symFun) subs(symFun,{X_dot_d,Y_dot_d,Z_dot_d,totD},...
        {'X_dot_d','Y_dot_d','Z_dot_d','totD'});
    
    centroidX=symfun(X_dot_d/totD,varArray);
    centroidY=symfun(Y_dot_d/totD,varArray);
    centroidZ=symfun(Z_dot_d/totD,varArray);
    
    centroid=symfun([centroidX,centroidY,centroidZ],varArray);
    lengthedgeJac=jacobian(centroid);
    lengthedgeHes=[hessian(centroidX),hessian(centroidY),hessian(centroidZ)];
    
    ccodeGen.f=custcccode(crep(centroid));
    ccodeGen.df=custcccode(crep(lengthedgeJac));
    if numel(ind)>2
        ccodeGen.ddf='';%custcccode(crep(lengthedgeHes));
    elseif numel(ind)==2
        mask=zeros(6);
        mask(1:3,4:6)=1;
        mask(4:6,1:3)=1;
        mask=repmat(mask,[1,3]);
        lengthedgeHes=lengthedgeHes.*mask;
        ccodeGen.ddf=custcccode(crep(lengthedgeHes));
    else
        ccodeGen.ddf=custcccode(crep(lengthedgeHes));
    end
    
    ccodeGen.symfun(1).fun=centroid;
    ccodeGen.symfun.jac=transpose(lengthedgeJac);
    ccodeGen.symfun.hes=crep(lengthedgeHes);
    
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
    
    ccodeGen=CCodeGenStruct('const vector<double>&','p1','p2','v0','v01','v02','v11','v12','c');
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
    
    ccodeGen.symfun(1).fun=lengthedge;
    ccodeGen.symfun.jac=lengthedgeJac;
    ccodeGen.symfun.hes=lengthedgeHes;
end

