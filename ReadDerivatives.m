function [coords, deriv]=ReadDerivatives(derivFile, numCoords)
    
    fid= fopen([derivFile,int2str(numCoords),'.csv'], 'r');
    coords = zeros(numCoords,3,0);
    deriv = zeros(3,3*numCoords,0);
    numDim = 1;
    numC = 1;
    numD = 1;
    while ~feof(fid)
        str = fgetl(fid);
        nums = str2num(str);
        if numel(nums)==0
            numDim=numDim+1;
            numC = 1;
            numD = 1;
        elseif numel(nums)==3
            coords(numC,:, numDim) = nums;
            numC=numC+1;
        else
            deriv(numD,:, numDim) = nums;
            numD=numD+1;
        end
    end
    fclose(fid);
    figure;
    ax=subplot(1,3,1);
    hold on;
    PlotDerivs(ax,numCoords,1,deriv);
    ax=subplot(1,3,2);
    hold on;
    PlotDerivs(ax,numCoords,2,deriv);
    ax=subplot(1,3,3);
    hold on;
    plot3(ax,coords([1:end,1],1,1),coords([1:end,1],2,1),coords([1:end,1],2,1))
    plot3(ax,coords([1:end,1],1,2),coords([1:end,1],2,2),coords([1:end,1],2,2))
    plot3(ax,coords([1:end,1],1,end),coords([1:end,1],2,end),coords([1:end,1],2,end))
    view(-35,20)
end

function []=PlotDerivs(ax,numCoords,step,deriv)
    nums = 1:numCoords:numCoords*3;
    for ii=1
        for jj=[nums,nums+step]
            plot(ax,squeeze(deriv(ii,jj,:)),'DisplayName',num2str([ii,jj]));
        end
    end
    legend()
end