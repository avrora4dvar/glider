%% Remove fields from files

removeList={'w','tke','AKv','AKt','AKs'}; 
expDir='/home/aruba/vol1/ipasmans/Exp/Exp51/'; 
expName='Exp51'; 


for expNo=2371:3:2392
    winDir=sprintf('%s/%s_%-.4d',expDir,expName,expNo);
    if ~exist(winDir,'dir'); continue; end
    display(winDir); 
    
    nlDir=sprintf('%s/%s_%-.4d/NL',expDir,expName,expNo);
    for iFile=1:99
        fname=sprintf('%s/ocean_his_%-.4d.nc',nlDir,iFile); 
        tmpname=sprintf('%s/tmp.nc',nlDir);
        if ~exist(fname,'file'); continue; end
        display(fname); 
        
        unix(sprintf('mv %s %s',fname,tmpname)); 
        reduce_his_file(tmpname,fname,removeList);
    end
    if exist(tmpname,'file'); delete(tmpname); end
    
    ensDir=sprintf('%s/%s_%-.4d/Ens',expDir,expName,expNo);
    for iFile=1:999
        fname=sprintf('%s/Member_%-.3d/ini_for.nc',ensDir,iFile); 
        tmpname=sprintf('%s/tmp.nc',expDir);
        if ~exist(fname,'file'); continue; end
        display(fname); 
        
        unix(sprintf('mv %s %s',fname,tmpname)); 
        reduce_his_file(tmpname,fname,removeList);
        
        fname=sprintf('%s/Member_%-.3d/ini.nc',ensDir,iFile); 
        tmpname=sprintf('%s/tmp.nc',expDir);
        if ~exist(fname,'file'); continue; end
        display(fname); 
        
        unix(sprintf('mv %s %s',fname,tmpname)); 
        reduce_his_file(tmpname,fname,removeList);
    end
    if exist(tmpname,'file'); delete(tmpname); end
    
end
