function []=MakeTexEquation(symeq,fileName,filePath)

    str=latex(symeq);

    str=regexprep(str,'\\_\\_','');
    str=regexprep(str,'\\_','_');
    sqrtCont=regexp(str,'\\sqrt{[^s]*}\^2}','match');
    sqrtCont=unique(sqrtCont);
    keepLog=cellfun(@numel,regexp(sqrtCont,'\{'))==...
        cellfun(@numel,regexp(sqrtCont,'\}'));
    sqrtCont=sqrtCont(keepLog);
    n=max(cellfun(@str2num,regexprep(unique(regexp(str,'_[0-9]*','match')),'_','')));
    if numel(sqrtCont)>=n
        sqrtCont=sqrtCont(1:n);
        sqrtCont=sqrtCont([1,3:end,2]);
    end
    translaTable=cell(0,3);

    ord=[1,n,2:n-1];
    fmatFun=@(dim,ind1,ind2) [dim,'_{',ind1,',',ind2,'}'];
    translaTable(end+1,1:2)={'X_dot_d',['\\sum x_{i,i+1}','d_i']};
    translaTable(end+1,1:2)={'Y_dot_d',['\\sum y_{i,i+1}','d_i']};
    translaTable(end+1,1:2)={'Z_dot_d',['\\sum z_{i,i+1}','d_i']};
    translaTable(end+1,1:2)={'totD',['\\sum d']};

    translaTable(end+1,1:2)={'\\left\(\\frac\{\\mathrm\{([xyz])_([0-9]*)\}\}\{2\} \+ \\frac\{\\mathrm\{[xyz]_([0-9]*)\}\}\{2\}\\right\)',fmatFun('$1','$2','$3')};
    translaTable(end,3)={false};
    for ii=1:numel(sqrtCont)
        d=['d_',int2str(ii)];
        %str=regexprep(str,regexptranslate('escape',sqrtCont{ii}),d);
        translaTable(end+1,1:2)={sqrtCont{ii},d};
    end
    for ii=1:numel(sqrtCont)
        d=['d_',int2str(ii),'^2'];
        translaTable(end+1,1:2)={sqrtCont{ii}(7:end-1),d};
    end
    % Simplify sums
    D='';
    for ii=ord
        d=['d_',int2str(ii)];
        if ii>1
            D=[D,' + '];
        end
        D=[D,d];
    end
    translaTable(end+1,1:2)={D,'\\sum d'};

    %str=regexprep(str,regexptranslate('escape',D),'\\sum d');
    for jj={'x','y','z'}
        D='';
        for ii=ord
            ind2=sort([ii,(mod(ii,n)+1)]);
            d=[fmatFun(jj{1},int2str(ind2(1)),int2str(ind2(2))),'\, ','d_',int2str(ii)];
            if ii>1
                D=[D,' + '];
            end
            D=[D,d];
        end
        %str=regexprep(str,regexptranslate('escape',D),['\\sum ',jj{1},'_{i+1}','d_i']);
        translaTable(end+1,1:2)={D,['\\sum ',jj{1},'_{i,i+1}','d_i']};
    end

    keepTrans=true([size(translaTable,1) 1]);
    for ii=1:size(translaTable,1)
        if isempty(translaTable{ii,3}) || translaTable{ii,3}
            str=regexprep(str,regexptranslate('escape',translaTable{ii,1}),translaTable{ii,2});
            if ii<size(translaTable,1)
                translaTable{ii,1}=[translaTable{ii,1},' \\'];
            end
            translaTable{ii,2}=[translaTable{ii,2},' = '];

        else
            str=regexprep(str,translaTable{ii,1},translaTable{ii,2});
            keepTrans(ii)=false;
        end
    end
    translaTable=translaTable(find(keepTrans),:);
    translaTable=[translaTable(:,2),translaTable(:,1)]';
    translaTable=regexprep(translaTable,'\\\\sum','\\sum');


    str=regexprep(str,'\&','\&\n');
    if exist('filePath','var');fileName=[filePath,filesep,fileName];end
    fid=fopen(fileName,'w');
    fprintf(fid, '%s',str);
    fclose(fid);
    fid=fopen(regexprep(fileName,'\.tex$','_trtble.tex'),'w');
    fprintf(fid, '%s',char(translaTable{:})');
    char(translaTable{:})
    fclose(fid);
end
