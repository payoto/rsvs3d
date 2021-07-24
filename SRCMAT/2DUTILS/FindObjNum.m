
function sub=FindObjNum(object,index,objInd)
    % The shape of the output is horribly inconsistant
    % finds the array index from a snaxel number
    if nargin==2
        objInd=[object(:).index];
    end

%     sub=FindObjNum_MexTest(index,objInd);
     %sub=FindObjNum_Matlab(index,objInd);

     % Mex
%     if isempty(objInd)
%         warning('list of object passed to findobjnum was empty')
%     end
    sub=FindObjNum_MEX(index,objInd,length(index),length(objInd));

    if numel(index)~=1
        sub=sub';
    end

  % Matlab
%     sub=zeros(length(index),1);
%     additionalSlots=0;
%     for ii=1:length(index)
%
%         snaxLog=objInd==index(ii);
%         jj=ii+additionalSlots;
%         subInter=find(snaxLog);
%         if isempty(subInter)
%             sub(jj)=0;
%         elseif numel(subInter)>1
%             sub(jj:jj+length(subInter)-1)=subInter;
%             additionalSlots=additionalSlots+numel(subInter)-1;
%         else
%             sub(jj)=subInter;
%
%         end
%     end


% Test match
%     if any(sub~=reshape(sub2,size(sub)))
%         error('subs are different')
%
%     end
end

function sub=FindObjNum_Matlab(index,objInd)

    sub=zeros(length(index),1);
    additionalSlots=0;
    for ii=1:length(index)

        snaxLog=objInd==index(ii);
        jj=ii+additionalSlots;
        subInter=find(snaxLog);
        if isempty(subInter)
            sub(jj)=0;
        elseif numel(subInter)>1
            sub(jj:jj+length(subInter)-1)=subInter;
            additionalSlots=additionalSlots+numel(subInter)-1;
        else
            sub(jj)=subInter;

        end
    end

end
function sub=FindObjNum_MexTest(index,objInd)

    sub=FindObjNum_MEX(index,objInd,length(index),length(objInd));

end
