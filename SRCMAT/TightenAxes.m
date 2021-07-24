function TightenAxes(ax, m, n, concatenateFigure)
    %moves the axes around to make total use of the available space.
    %ax is an array of axes handles, which assumes that we count accross left to
    %m is the number of rows
    %n is the number of columns
    %Each axes and it's labels is allowed to have 1/n width and 1/m height which
    %drives the width and height of the axes themselves.
    %right and then down. i.e. m = 2, n = 3
    %	ax(1) - left top
    %	ax(2) - middle top
    %	ax(3) - right top
    %	ax(4) - left bottom
    %	etc
    % copyright Matt Edwards 2019

    if ~exist('concatenateFigure','var')
        concatenateFigure =false;
    end

    mrgn = 0.005; %effectively applied to the right and to the top of each axes
    ax = findobj(ax, 'type', 'axes');
    %check that the dimensions and number of axes is equivalent
    assert(~isempty(ax), 'No axes were found at the targets provided.')
    if ischar(m)
        [m,n] = CalculateLayout(m, numel(ax));
    end
    if concatenateFigure
        [ax] = ConcatenateFigures(ax);
    end
    %initialise w and h to their largest possible values
    w = 1;
    h = 1;
    hOffset = 0;
    vOffset = 0;
    %work out what width and height we are allowed
    for ii = 1:numel(ax)
        %each axes is allowed 1/n of the width. Therefore we need to work out what
        %the width needs to be such that they are all the same, when they are allowed
        %to be padded by their TightInset (i.e. the axis ticks)
        w_bf = (1/n) - ax(ii).TightInset(1) - ax(ii).TightInset(3) - mrgn;
        if w_bf < w
            w = w_bf;
        end
        h_bf = (1/m) - ax(ii).TightInset(2) - ax(ii).TightInset(4) - mrgn;
        if h_bf < h
            h = h_bf;
        end
        if ax(ii).TightInset(1) > hOffset
            hOffset = ax(ii).TightInset(1);
        end
        if ax(ii).TightInset(2) > vOffset
            vOffset = ax(ii).TightInset(2);
        end
    end

    %now loop over the axes calculating their nn and mm index position in
    % the matrix of figures and aligning the axes.
    for ii = 1:numel(ax)
        nn = mod(ii-1,n)+1;
        mm = 1+floor((ii-1)/n);
        left = (nn-1)*(1/n)+hOffset;
        bottom = (m-mm)*(1/m) + vOffset;
        ax(ii).Position = [left, bottom, w, h];
    end

end

function [ax] = ConcatenateFigures(ax)

    sameParents = true;
    ii = 1;
    while sameParents && ii <= numel(ax)
        sameParents = ax(ii).Parent == ax(1).Parent;
        ii=ii+1;
    end
    if ~sameParents
        parents = unique([ax.Parent]);
        newFigName = parents(1).Name;
        sizeFig = parents(1).Position(3:4);
        for ii = 2:numel(parents)
            newFigName = [' + ', parents(ii).Name];
            sizeFig = max(parents(ii).Position(3:4),sizeFig);
        end
        fig = figure('Name',newFigName, 'position', [200,200, sizeFig]);
        ax = copyobj(ax, fig);
    end

end

function [m,n] = CalculateLayout(charIn, tot)
    % automatically calculates the layout

    switch charIn
        case {'auto', 'rowmajor'}
            n = ceil(sqrt(tot));
            m = ceil(tot/n);
        case 'colmajor'
            m = ceil(sqrt(tot));
            n = ceil(tot/m);
        otherwise
            error('Unknown calculated layout.')
    end
end
