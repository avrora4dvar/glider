%relabel ocean_his_ ocean_avg files

roms_dir='/home/aruba/vol2/ipasmans/Exp/Exp26g_20110425'; 
dateRef_roms=datenum('2005-01-01'); 
dateRef_dir=datenum('2005-01-01'); 
%% find files

con=dir(roms_dir); 

for iF=1:length(con)
   if strncmpi(con(iF).name,'ocean_his',9)
       fileIn=fullfile(roms_dir,con(iF).name); 
       t=ncread(fileIn,'ocean_time',1,1); 
       
       t=t/24/3600+dateRef_roms; 
       t=floor(t-dateRef_dir); 
       
       fileOut=fullfile(roms_dir,sprintf('ocean_his_%-.4d.nc',t));
       command=sprintf('mv %s %s',fileIn,fileOut); 
       unix(command); 
   end
end

for iF=1:length(con)
   if strncmpi(con(iF).name,'ocean_avg',9)
       fileIn=fullfile(roms_dir,con(iF).name); 
       t=ncread(fileIn,'ocean_time',1,1); 
       
       t=t/24/3600+dateRef_roms; 
       t=floor(t-dateRef_dir); 
       
       fileOut=fullfile(roms_dir,sprintf('ocean_avg_%-.4d.nc',t));
       command=sprintf('mv %s %s',fileIn,fileOut); 
       unix(command); 
   end
end
