%% remove uncessary files from exp dirs

romsDir='/home/aruba/vol2/ipasmans/Exp/Exp35/';
expName='Exp36'; 
winNo=[2392:3:2413]; 

%current directory
nowDir=cd(); 

for iWin=winNo
   %directory to be cleaned
   expDir=fullfile(romsDir,sprintf('%s_%-.4d',expName,iWin)); 
   if ~exist(expDir,'dir'); continue; end; 
   display(expDir); 
   
   %clean Iter0
   fNames={'back.nc'}; 
   iterDir=fullfile(expDir,'Iter0'); 
   if exist(iterDir,'dir')
       cd(iterDir)
       for k=1:length(fNames)
           if exist(fNames{k},'file'); delete(fNames{k}); end
       end
   end
   
   %clean Iter1
   fNames={'r.nc','per_ana.nc'};
   fNames={}; 
   iterDir=fullfile(expDir,'Iter1'); 
   if exist(iterDir,'dir')
       cd(iterDir)
       for k=1:length(fNames)
           if exist(fNames{k},'file'); delete(fNames{k}); end
       end
   end
   
   %clean nonlinear
   fNames={'ini.nc'}; 
   iterDir=fullfile(expDir,'NL'); 
   if exist(iterDir,'dir')
       cd(iterDir)
       for k=1:length(fNames)
           if exist(fNames{k},'file'); delete(fNames{k}); end
       end
   end

   
   %clean ensemble
   ensDir=fullfile(expDir,'Ens'); 
   if exist(ensDir,'dir')
       
       fNames={}; 
       cd(ensDir)
       for k=1:length(fNames)
           if exist(fNames{k},'file'); delete(fNames{k}); end
       end
       
       fNames={'ocean_avg_for.nc','ocean_his_for.nc','per_ana.nc','Ad.nc'}; 
       for iMember=1:999
           memberDir=fullfile(ensDir,sprintf('Member_%-.3d',iMember));
           if ~exist(memberDir,'dir'); break; end
           cd(memberDir); 
           for k=1:length(fNames)
               if exist(fNames{k},'file'); delete(fNames{k}); end
           end 
       end
       
   end
   
end

%back to current dir
cd(nowDir); 