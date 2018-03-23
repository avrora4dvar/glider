%% Calculate array modes
clear all; clc; 
addpath('/home/server/student/homes/ipasmans/Matlab/figures/obslist'); 

binDir='/home/aruba/vol2/ipasmans/Bin_Exp35'; 
romsDir='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_2395'; 

%% Find ritz vectors

r=ncread(fullfile(romsDir,'Iter1','pcg.nc'),'r');
Br=ncread(fullfile(romsDir,'Iter1','pcg.nc'),'Br'); 


[pcg,iter,obs]=iter_data(fullfile(romsDir,'Iter1','pcg.nc')); 

%% Construct the Ritz vectors

ritzV=iter{end}.lanczos_vector;
ritzS=iter{end}.ritz;

r=r(:,1:length(ritzS)); 
Br=Br(:,1:length(ritzS)); 

for k=1:size(ritzV,2)
    eigV(:,k)=r*ritzV(:,k); 
    eigBV(:,k)=Br*ritzV(:,k); 
    eigNorm(k)=sqrt(eigV(:,k)'*eigBV(:,k)); 
    eigBV(:,k)=eigBV(:,k).*obs.sig/eigNorm(k); 
    eigV(:,k)=eigV(:,k)/eigNorm(k); 
end

nVec=find(cumsum(ritzS)/sum(ritzS)<.9,1,'last');
nVec=56;

%% Convert from dual to primal

outDir=fullfile(romsDir,'Ritz'); 
%if exist(outDir,'dir'); unix(sprintf('rm -r %s',outDir)); end
%mkdir(outDir); 

%Background
if ~exist(fullfile(romsDir,'Iter0','back.nc'),'file')
    fileIn=fullfile(romsDir,'Iter0','back.in');
    fileOut=fullfile(romsDir,'Iter0','back.out');
    if exist(fileOut,'file'); delete(fileOut); end
    unix(sprintf('%s/create_background < %s > %s',...
        binDir,fileIn,fileOut));
end

for k=47:nVec
   outDir=fullfile(romsDir,'Ritz',sprintf('Member_%-.3d',k)); 
   if ~exist(outDir,'dir'); mkdir(outDir); end
   
   fname=fullfile(outDir,'v.nc'); 
   copyfile(fullfile(romsDir,'Iter0','sample_for.nc'),fname); 
   ncwrite(fname,'K',eigV(:,k)./obs.sig);
   pause(1);
   
   %Ad_avrora
   display(sprintf('Starting ad_avrora %d',k)); 
   copyfile(fullfile(romsDir,'case_def.txt'),fullfile(outDir,'case_def.txt')); 
   fid=fopen(fullfile(outDir,'case_def.txt'),'a');
   fprintf(fid,'tiles=(%d %d)\n',2,3); 
   fprintf(fid,'background_file=%s\n',fullfile(romsDir,'Iter0','back.nc')); 
   fprintf(fid,'input=%s\n',fname); 
   fprintf(fid,'output=%s\n',fullfile(outDir,'Adv.nc')); 
   fclose(fid); 
   unix(sprintf('sh %s/run_ad.sh %s %s %s &',binDir,...
       fullfile(outDir,'case_def.txt'),...
       outDir,'aruba2')); 
   
   %Wait
   while true
       if mod(k-46,4)~=0 && k<nVec; break; end
       flag_go=1;   
       for l=k:-1:1
           fname=fullfile(romsDir,'Ritz',sprintf('Member_%-.3d',l),'ad_interp.out');
           if ~exist(fname,'file')
               flag_go=0;
           else
               [s,cmdout]=unix(sprintf('tail -n 1 %s',fname));
               flag_go=flag_go*~isempty(regexp(cmdout,'DONE','once'));
           end
       end
       pause(5);
       if flag_go~=0; break; end
   end
      
   
end

for k=47:nVec
    outDir=fullfile(romsDir,'Ritz',sprintf('Member_%-.3d',k)); 
    
   %Covariance
   display(sprintf('Starting covariance %d',k)); 
   copyfile(fullfile(romsDir,'case_def.txt'),fullfile(outDir,'case_def.txt')); 
   fid=fopen(fullfile(outDir,'case_def.txt'),'a');
   fprintf(fid,'background_file=%s\n',fullfile(romsDir,'Iter0','back.nc')); 
   fprintf(fid,'input=%s\n',fullfile(outDir,'Adv.nc')); 
   fprintf(fid,'output=%s\n',fullfile(outDir,'Cv.nc')); 
   fclose(fid); 
   unix(sprintf('sh %s/run_cov_bal.sh %s %s %s &',binDir,...
       fullfile(outDir,'case_def.txt'),...
       outDir,'aruba2')); 
   
   %Wait 
   while true
       if mod(k-46,4)~=0 && k<nVec; break; end
       flag_go=1;  
       for l=k:-1:1
           fname=fullfile(romsDir,'Ritz',sprintf('Member_%-.3d',l),'cov_tl.out');
           if ~exist(fname,'file')
               flag_go=0;
           else
               [s,cmdout]=unix(sprintf('tail -n 1 %s',fname));
               flag_go=flag_go*~isempty(regexp(cmdout,'output','once'));
           end
       end
       pause(5);
       if flag_go~=0; break; end
   end
   
end

for k=47:nVec
    outDir=fullfile(romsDir,'Ritz',sprintf('Member_%-.3d',k)); 
    delete(fullfile(outDir,'Adv.nc')); 
end