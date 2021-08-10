function [] = include_Validation()
    %FUNCTIONLIST allows local functions to be used globally once it has
    %been used.

    funcHandles=localfunctions;
    funcDir=[cd,'\Automated_Function_Directory'];
    HeaderActivation(funcHandles,funcDir)

end

%%

function [outLine]=GenerateValidationData(procdat)

    kk=1;
    outLine{kk,1}=sum([procdat(:).warningnum]);outLine{kk,2}='int';outLine{kk,3}='warnings';
    kk=kk+1;
    outLine{kk,1}=sum([procdat(:).errornum]);outLine{kk,2}='int';outLine{kk,3}='errors';
    kk=kk+1;
    outLine{kk,1}=[sum([procdat(:).termination]),length(procdat)];
    outLine{kk,2}='int';outLine{kk,3}='termination';
    kk=kk+1;
    outLine{kk,1}=sum([procdat(:).cputime]);outLine{kk,2}='double';outLine{kk,3}='total cputime';
    kk=kk+1;
    outLine{kk,1}=sum([procdat(:).snaketime]);outLine{kk,2}='double';outLine{kk,3}='total snaketime';
    kk=kk+1;
    outLine{kk,1}=mean([procdat(:).volerror]);outLine{kk,2}='double';outLine{kk,3}='vol err - Mean';
    kk=kk+1;
    outLine{kk,1}=std([procdat(:).volerror]);outLine{kk,2}='double';outLine{kk,3}='vol err - STD';
    kk=kk+1;
    outLine{kk,1}=max([procdat(:).volerror]);outLine{kk,2}='double';outLine{kk,3}='vol err - max';
    kk=kk+1;
    outLine{kk,1}=min([procdat(:).volerror]);outLine{kk,2}='double';outLine{kk,3}='vol err - min';
    kk=kk+1;
    outLine{kk,1}=mean([procdat(:).velerror]);outLine{kk,2}='double';outLine{kk,3}='vel err - Mean';
    kk=kk+1;
    outLine{kk,1}=std([procdat(:).velerror]);outLine{kk,2}='double';outLine{kk,3}='vel err - STD';
    kk=kk+1;
    outLine{kk,1}=max([procdat(:).velerror]);outLine{kk,2}='double';outLine{kk,3}='vel err - max';
    kk=kk+1;
    outLine{kk,1}=min([procdat(:).velerror]);outLine{kk,2}='double';outLine{kk,3}='vel err - min';
    kk=kk+1;
    outLine{kk,1}=sum([procdat(:).length]);outLine{kk,2}='double';outLine{kk,3}='total length';

end


function [sumLines]=GenerateSummaryLine(procdat)

    [outLine]=GenerateValidationData(procdat);

    separator='';
    sumLines{1}='# ';
    sumLines{2}='';
    for ii=1:size(outLine,1)
        sumLines{1}=[sumLines{1},separator,outLine{ii,3}];
        sumLines{2}=[sumLines{2},separator,ProcesstoString(outLine{ii,1},...
            outLine{ii,2})];
        separator=' , ';
    end


end

function [procdat]=ProcDatTemplate()

    procdat=struct('caseName','','warningnum',0,'errornum',0,'cputime',0,'termination',false, 'niter',0,...
            'snaketime',0,'volerror',1,'velerror',1,'length',0,'path','');

end


function [strOut]=ProcesstoString(inputVar,classin)

    if nargin==1
        classin=class(inputVar);
    end

    switch classin
        case 'int'
            strOut=int2str(inputVar);
        case 'logical'
            strOut=int2str(inputVar);
        case 'double'
            if all(mod(inputVar,1)==0)
                strOut=int2str(inputVar);
            else
                strOut=num2str(inputVar,' %30.24e ');
            end
        case 'char'
            strOut=inputVar;
        otherwise
            strOut='';

    end

end
