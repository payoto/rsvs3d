%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              University of Bristol
%      Department of Aerospace Engineering
%                     2015
%
%          Subdivision of Surfaces
%      for Aerodynamic shape parametrisation
%          - WRAPPER for Mex Code -
%             Alexandre Payot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = include_Mex_Wrapper()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.

    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)

end

%% Mex Wrapper functions

function []=GridInit_Wrapper()


    GridInit_MEX

    clear GridInit_MEX

end

function  [refinedGrid,connec]=GridRefine_Wrapper(grid,cellRefInd,refineSplit)

    % Minimize work by modifying position vector in grid?
    refineVecs=vertcat(grid.cell(:).refineVec);

    [grid.cell(:).refineVec]=deal(grid.cell(:).index);
    [grid.cell(:).refineLvl]=deal(1);
    sizCell=size(cellRefInd);
    nCell=length(grid.cell);
    nEdge=length(grid.edge);
    nVert=length(grid.vertex);
    nLevels=2;
    levelSize=[1 nCell refineSplit];

    if sizCell(2)==1
        cellRefPos=FindObjNum([],cellRefInd,[grid.cell(:).index])-1;
    elseif sizCell(2)==2
        cellRefPos=cellRefInd(:,2)-1;
        cellRefInd=cellRefInd(:,1);
    else
        error('Cell Not right');

    end
    if sum(cellRefPos<0)
        error('Invalid Cell Specified for refinement');
    end
    nRefine=length(cellRefInd);
    [refinedGrid,connec]=GridRefine_MEX(nLevels,levelSize,...
                    nCell,nEdge,nVert,grid,...
                    nRefine,cellRefInd,cellRefPos);
%     GridRefine_MEX(nLevels,levelSize,...
%                     nCell,nEdge,nVert,grid,...
%                     nRefine,cellRefInd,cellRefPos);
%                 refinedGrid=grid;
%                 connec=1;
    clear GridRefine_MEX
    if ~exist('refinedGrid','var'),error('refinedGrid not assigned during call to GridREfine_MEX');end
    if ~exist('connec','var'),conec=[];warning('connec not assigned during call to GridREfine_MEX');end
    %refinedGrid=grid;

end
