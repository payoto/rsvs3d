function [spaceMat]=VenusEarthOrbitImage(n,p1,r1,p2,r2,nRot)
    
    spaceMat=zeros(2*n+1);
    isInSpace=zeros(2*n+1);
    [indMat1,indMat2]=meshgrid(1:2*n+1,1:2*n+1);
    subList=[indMat1(:),indMat2(:)]';
    rmax=max(r1,r2);
    r1=r1/rmax;
    r2=r2/rmax;
    tperRot=50;
    dt=min(p1,p2)/tperRot;
    cIncr=0.1/1.8/10.4;
    t=0;
    plotPlanets=@(p1,p2) plot([p1(1),p2(1)],[p1(2),p2(2)],'wo-','markersize',3,'MarkerFaceColor',[1 1 1]);
    % spaceMat(n+1,n+1)=(0,0)
    % spaceMat(2n+1,2n+1)=(1,1)
    % spaceMat(1,1)=(-1,-1)
    kk=1;
    pos1=[r1;0];
    pos2=[r2;0];
    h=figure;
    while t<nRot*p1
%         cBound=[min(spaceMat(spaceMat(:)>0)) max(spaceMat(:))];
%         if numel(cBound)==1
%             cBound=[0 1];
%         end
        imshow(spaceMat,[])
        
%         plotPlanets(pos1(:,kk),pos2(:,kk));
%         hold off
        movFrame(kk)=getframe(h);
        tSave(kk)=t;
        
        kk=kk+1;
        %image(spaceMat)
        [pos1(:,kk)]=PositionOrbit(pos1(:,kk-1),p1,dt);
        [pos2(:,kk)]=PositionOrbit(pos2(:,kk-1),p2,dt);
        [poly]=FindPoly(pos1(:,kk-1:kk),pos2(:,kk-1:kk));
        [spaceMat]=LogicalInPoly(poly,isInSpace,subList,n)*cIncr+spaceMat;
        %plot(poly(1,:),poly(2,:),'k-')
        t=t+dt;
        tSave(kk)=t;
    end
        imshow(spaceMat,[])
        
%         plotPlanets(pos1(:,kk),pos2(:,kk));
%         hold off
        movFrame(kk)=getframe(h);
%     figure, hold on
%     plot(pos1(1,:),pos1(2,:))
%     plot(pos2(1,:),pos2(2,:))
%     figure, hold on
%     plot(tSave,pos1(1,:),tSave,pos1(2,:))
%     plot(tSave,pos2(1,:),tSave,pos2(2,:))
    if 1
        fps=30;
        quality=100;
        MakeVideo(movFrame,fps,quality,'planets.avi');
    end
end

function [newpos]=PositionOrbit(pos,p,dt)
    
    rot=[cos(dt/p*2*pi) -sin(dt/p*2*pi);
        sin(dt/p*2*pi) cos(dt/p*2*pi)];
    newpos=rot*pos;
    
end

function [poly]=FindPoly(pos1,pos2)
    rot90=[0 -1;1 0];
    
%      if dot((pos1(:,2)-pos1(:,1))./sqrt(sum((pos1(:,2)-pos1(:,1)).^2)),...
%              (pos2(:,2)-pos2(:,1))./sqrt(sum((pos2(:,2)-pos2(:,1)).^2)))<0
% 
%     
% 
%         mat=[(rot90*(pos2(:,2)-pos1(:,2))),(rot90*(pos1(:,1)-pos2(:,1)))]';
%         
%         vec(1,1)=(rot90*(pos2(:,2)-pos1(:,2)))'*pos1(:,2);
%         vec(2,1)=(rot90*(pos1(:,1)-pos2(:,1)))'*pos2(:,1);
%         
%         xIntersect=mat\vec;
%         poly{1}=([pos1,xIntersect]);
%         poly{2}=([pos2,xIntersect]);
%         
%     else
    if xor(dot((pos2(:,2)-pos2(:,1)),rot90*(pos2(:,2)-pos1(:,2)))<0,...
            dot((pos1(:,2)-pos1(:,1)),rot90*(pos2(:,2)-pos1(:,2)))<0) || ...
            xor(dot((pos2(:,1)-pos2(:,2)),rot90*(pos2(:,1)-pos1(:,1)))<0,...
            dot((pos1(:,1)-pos1(:,2)),rot90*(pos2(:,1)-pos1(:,1)))<0)
        mat=[(rot90*(pos2(:,2)-pos1(:,2))),(rot90*(pos1(:,1)-pos2(:,1)))]';
        
        vec(1,1)=(rot90*(pos2(:,2)-pos1(:,2)))'*pos1(:,2);
        vec(2,1)=(rot90*(pos1(:,1)-pos2(:,1)))'*pos2(:,1);
        
        xIntersect=mat\vec;
        poly{1}=([pos1,xIntersect]);
        poly{2}=([pos2,xIntersect]);
     else
        
        kk=2;
        poly{1}=[pos1(:,kk-1:kk),pos2(:,kk:-1:kk-1)];
        %poly{2}=flip([pos1(:,kk-1:kk),pos2(:,kk:-1:kk-1)]);
        %poly{2}=[pos2(:,kk-1:kk),pos1(:,kk:-1:kk-1)];
    end
    
end

function [isInSpace]=LogicalInPoly(poly,isInSpace,subList,n)
    
    rot90=-[0 -1;1 0];
    
    logCond=false([numel(poly),size(subList,2)]);
    for ii=1:numel(poly)
        
        normVecs=rot90*(poly{ii}(:,[2:end,1])-poly{ii});
        condMat=normVecs'*(subList-n+1)*2/(2*n+1);
        conds=sum(normVecs.*poly{ii},1);
        for jj=1:numel(conds)
            condMat(jj,:)=condMat(jj,:)<=conds(jj);
        end
        logCond(ii,:)=all(condMat,1) | ~any(condMat,1);
        
    end
    if sum(sum(logCond))>0
        if numel(poly)>1
            %disp('multipoly works')
        end
    else
        if numel(poly)>1
            disp('multipoly not working')
        else
            plotPoints= @(points) plot(points([1:end],1),points([1:end],2),'*-');
            plotPoints(poly{1}')
            disp('single not working')
        end
    end
    isInSpace(any(logCond,1))=1;
end