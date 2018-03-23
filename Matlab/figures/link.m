% create soft links to analysis and forecasts
clear all; clc; 

obsDir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp37'; 
romsDir='/home/aruba/vol2/ipasmans/Exp/Exp35/';
expName='Exp42'; 
forDir='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp42_for/'; 
anaDir='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp42_ana/';

tl_sample='/home/aruba/vol2/ipasmans/Bin_Exp35/tl_sample'; 
grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
dateRef=datenum('2005-01-01'); 
winLength=3; 
flag_ens=false; 
flag_sample=true;

%% create output directories

if ~exist(forDir,'dir')
    mkdir(forDir);
end
if ~exist(anaDir,'dir')
    mkdir(anaDir); 
end

%% make links

flag_for1=true; 
con=dir(romsDir);

for iDir=1:length(con)
   winNo=regexp(con(iDir).name,[expName,'_(\d+)'],'tokens'); 
   
   if length(winNo)~=1
       continue; 
   else
       winNo=str2num(winNo{1}{1}); 
   end
   
   nlDir=fullfile(romsDir,con(iDir).name,'NL'); 
   if ~exist(nlDir,'dir'); continue; end
   
   conNl=dir(nlDir); 
   for iFile=1:length(conNl)
       hisNo=regexp(conNl(iFile).name,'ocean_his_(\d+).nc','tokens'); 
       if isempty(hisNo); continue; end
       filename=fullfile(nlDir,conNl(iFile).name); 
       
        t=ncread(filename,'ocean_time',1,1);
        t=t/3600/24+dateRef;
        t=floor(t-dateRef); 
        
        
        if t>=winNo && t<winNo+winLength
            display(sprintf('linking %s to analysis',filename)); 
            %analysis
            command=sprintf('ln -snf -T %s %s',...
                filename,...
                fullfile(anaDir,sprintf('ocean_his_%-.4d.nc',t))); 
            unix(command); 
        else
            display(sprintf('linking %s to forecast',filename)); 
            %forecast
            command=sprintf('ln -snf -T %s %s',...
                filename,...
                fullfile(forDir,sprintf('ocean_his_%-.4d.nc',t))); 
            unix(command); 
        end
        
        if flag_for1 && t==winNo+winLength-1
            display(sprintf('linking %s to forecast',filename)); 
            flag_for1=false;
            %forecast
            command=sprintf('ln -snf -T %s %s',...
                filename,...
                fullfile(forDir,sprintf('ocean_his_%-.4d.nc',t))); 
            unix(command); 
        end
            
   end
   
   %% average
   
   conNl=dir(nlDir);
   for iFile=1:length(conNl)
       hisNo=regexp(conNl(iFile).name,'ocean_avg_(\d+).nc','tokens');
       if isempty(hisNo); continue; end
       filename=fullfile(nlDir,conNl(iFile).name);
       
       t=ncread(filename,'ocean_time',1,1);
       t=t/3600/24+dateRef;
       t=floor(t-dateRef);
       
       
       if t>=winNo && t<winNo+winLength
           display(sprintf('linking %s to analysis',filename));
           %analysis
           command=sprintf('ln -snf -T %s %s',...
               filename,...
               fullfile(anaDir,sprintf('ocean_avg_%-.4d.nc',t)));
           unix(command);
       else
           display(sprintf('linking %s to forecast',filename));
           %forecast
           command=sprintf('ln -snf -T %s %s',...
               filename,...
               fullfile(forDir,sprintf('ocean_avg_%-.4d.nc',t)));
           unix(command);
       end
   end
   
   %% Ensemble
   
   if flag_ens
       for iMem=1:39
           memDir=fullfile(romsDir,con(iDir).name,'Ens',sprintf('Member_%-.3d',iMem));
           memForDir=fullfile(forDir,sprintf('Member_%-.3d',iMem)); 
           memAnaDir=fullfile(anaDir,sprintf('Member_%-.3d',iMem)); 
           
           if ~exist(memForDir,'dir'); mkdir(memForDir); end
           if ~exist(memAnaDir,'dir'); mkdir(memAnaDir); end
           
           display(sprintf('Linking %s',memDir)); 
           command=sprintf('ln -snf -T %s %s',...
               fullfile(memDir,'ini_for.nc'),fullfile(memForDir,sprintf('ocean_his_%-.4d.nc',winNo))); 
           unix(command); 
           command=sprintf('ln -snf -T %s %s',...
               fullfile(memDir,'ini.nc'),fullfile(memAnaDir,sprintf('ocean_his_%-.4d.nc',winNo))); 
           unix(command); 
               
       end
   end
   
   %% Sample
   
   if flag_sample
    
       
       %sample ana
       display('Sampling analysis'); 
       fid=fopen(fullfile(nlDir,'sample.in'),'w+');
       fprintf(fid,'#\n%s\n',fullfile(obsDir,sprintf('obslist%-.4d.nc',winNo)));
       fprintf(fid,'#\n%s\n%s\n',nlDir,'ocean_his');
       fprintf(fid,'#\n%s\n',grdFile);
       fprintf(fid,'#\n%s\n%s\n',nlDir,'ocean_his');
       fprintf(fid,'#\n%s\n',fullfile(anaDir,sprintf('sample%-.4d.nc',winNo)));
       fprintf(fid,'#\n%s\n','K');
       fprintf(fid,'#\n%d,%d,%d',winNo,winNo+3,-winNo);
       fclose(fid);
       
       delete(fullfile(nlDir,'sample.out'));
       command=sprintf('%s < %s > %s',tl_sample,...
           fullfile(nlDir,'sample.in'),...
           fullfile(nlDir,'sample.out'));
       unix(command);
       
       %sample for
       display('Sampling forecast'); 
       fid=fopen(fullfile(nlDir,'sample.in'),'w+');
       fprintf(fid,'#\n%s\n',fullfile(obsDir,sprintf('obslist%-.4d.nc',winNo+3)));
       fprintf(fid,'#\n%s\n%s\n',nlDir,'ocean_his');
       fprintf(fid,'#\n%s\n',grdFile);
       fprintf(fid,'#\n%s\n%s\n',nlDir,'ocean_his');
       fprintf(fid,'#\n%s\n',fullfile(forDir,sprintf('sample%-.4d.nc',winNo+3)));
       fprintf(fid,'#\n%s\n','K');
       fprintf(fid,'#\n%d,%d,%d',winNo+3,winNo+6,-winNo-3);
       fclose(fid);
       
       delete(fullfile(nlDir,'sample.out'));
       command=sprintf('%s < %s > %s',tl_sample,...
           fullfile(nlDir,'sample.in'),...
           fullfile(nlDir,'sample.out'));
       unix(command);
       
   else
       %analysis
       command=sprintf('ln -snf -T %s %s',...
           fullfile(romsDir,con(iDir).name,'Iter1','sample_ana.nc'),...
           fullfile(anaDir,sprintf('sample%-.4d.nc',winNo))); 
       unix(command);
       %Forecast
       command=sprintf('ln -snf -T %s %s',...
           fullfile(romsDir,con(iDir).name,'Iter0','sample_for.nc'),...
           fullfile(forDir,sprintf('sample%-.4d.nc',winNo))); 
       unix(command);           
   end
   
end