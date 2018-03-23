function [info] = reduce_his_file(fileIn,fileOut,removeList)
%REDUCE_HIS_FILE fileOut=reduce_his_file(fileIn)
%   

%Read writeschema 
info=ncinfo(fileIn); 

%Possible fields to remove: w, AKv, AKt, AKs, tke, shflux,

%Copy only variables not in removeList
varList=[];
for iVar=1:length(info.Variables)
    if ~ismember(info.Variables(iVar).Name,removeList)
        varList=[varList,info.Variables(iVar)]; 
    end
end
info.Variables=varList; 

%Create new file
ncwriteschema(fileOut,info); 

%Copy fields
for iVar=1:length(info.Variables)
    val=ncread(fileIn,info.Variables(iVar).Name); 
    ncwrite(fileOut,info.Variables(iVar).Name,val); 
end

end

