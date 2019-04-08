function [] = ASOV3_interface()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.
    
    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)
    
end


%% ASO V3 interface

% Generate geom object
function [geom]=Grid2Geom3(grid)
    
    vertInd= [grid.vert.index];
    edgeInd= [grid.edge.index];

    faceInput = cell([numel(grid.surf),1]);
    for ii = 1: numel(grid.surf)
        [currVertSub]=OrderSurfaceVertList(grid, ii, edgeInd, vertInd);

        faceInput{ii} = currVertSub([1:end])';
        if grid.surf(ii).voluind(2)==0
            faceInput{ii} = flip(faceInput{ii});
        end
    end
    
    vertCoord = vertcat(grid.vert.coord);
    geom = geom3.UnstrucSurf(faceInput, vertCoord);
end

% 