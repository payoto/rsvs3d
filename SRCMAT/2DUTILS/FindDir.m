function [returnPath,returnName]=FindDir(rootDir,strDir,isTargDir)
    returnPath={};
    returnName={};
%     if iscell(rootDir)
%         subDir=dir(rootDir{1});
%         subDir(1:2)=[];
%         for ii=2:numel(rootDir)
%             partsubDir=dir(rootDir{ii});
%             partsubDir(1:2)=[];
%             subDir=[subDir,partsubDir];
%         end
%     else
        subDir=dir(rootDir);
        subDir(1:2)=[];
%     end
    nameCell={subDir(:).name};
    isprofileDirCell=strfind(nameCell,strDir);
    for ii=1:length(subDir)
        subDir(ii).isProfile=(~isempty(isprofileDirCell{ii})) && ...
            ~xor(subDir(ii).isdir,isTargDir);
    end
    
    returnSub=find([subDir(:).isProfile]);
    
    
    if isempty(returnSub)
        fprintf('FindDir Could not find requested item %s in:\n%s \n',strDir,rootDir)
    end
    for ii=1:length(returnSub)
        returnPath{ii}=[rootDir,filesep,subDir(returnSub(ii)).name];
        returnName{ii}=subDir(returnSub(ii)).name;
        
    end
    
    
    
end

