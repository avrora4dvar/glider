%% Get cross from model 
clear all; clc; 

addpath('..'); 
addpath('/home/server/student/homes/ipasmans/Matlab/roms/seawater/'); 
addpath('/home/server/student/homes/ipasmans/Matlab/figures/io/'); 

%directory with glider files
obsDir='/home/server/student/homes/ipasmans/Obs/glider/osu/';

%grid file
grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';


romsDir={'/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp40_ana'}; 


%time boundaries
dateRef=datenum('2005-01-01'); 
tLim=[2392,2413]+dateRef; 
dLim=[0,250]; 


%output
outFile='/home/server/student/homes/ipasmans/Data/article_glider/TS_alongglider_ana.mat';
clear model;

%% Grid

grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.z=ncread(grdFile,'z0_r'); 
grd.h=ncread(grdFile,'h'); 

%% read observations

obs=io_glider_osu(obsDir,'dir'); 
obs.source='http://gliderfs2.coas.oregonstate.edu/gliderweb/archive/gridded/slocum_gridded_data.php'; 
in=obs.t>=tLim(1) & obs.t<=tLim(2); 
obs.lon=obs.lon(in); obs.lat=obs.lat(in); obs.t=obs.t(in); 
obs.temp=obs.temp(in,:); obs.salt=obs.salt(in,:); obs.p=obs.p(in,:); obs.depth=obs.depth(in,:); 
obs.sig_temp=obs.sig_temp(in,:); obs.sig_salt=obs.sig_salt(in,:); obs.sig_p=obs.sig_p(in,:); 


%% Read model

model=[]; 
for iMod=1:length(romsDir)
   con=dir(romsDir{iMod}); 
   
   model1=struct('lon',[],'lat',[],'t',[],'depth',[],'temp',[],'salt',[],'dir',[]); 
   for iFile=1:length(con)
      if isempty(regexp(con(iFile).name,'ocean_avg_(\d+).nc','once')); continue; end
      fname=fullfile(romsDir{iMod},con(iFile).name); 
      display(fname); 
      
      tFile=ncread(fname,'ocean_time')/24/3600+dateRef;
      datestr(tFile)
      tLim=[floor(min(tFile)),ceil(max(tFile))]; 
      in=obs.t>tLim(1) & obs.t<tLim(2);
      if ~any(in); continue; end
      
      lat1=mean(obs.lat(in)); 
      [~,ilat1]=min( abs(grd.lat(1,:)-lat1) ); 
      lat1=grd.lat(1,ilat1); 
      
      temp1=ncread(fname,'temp',[1,ilat1,1,1],[Inf,1,Inf,Inf]); 
      salt1=ncread(fname,'salt',[1,ilat1,1,1],[Inf,1,Inf,Inf]); 
      temp1=mean(temp1,4); 
      salt1=mean(salt1,4); 
      
      model1.lon=[model1.lon;reshape(grd.lon(:,1),[],1)]; 
      model1.lat=[model1.lat;ones(size(grd.lon,1),1)*lat1]; 
      model1.t=[model1.t;ones(size(grd.lon,1),1)*mean(tLim)]; 
      model1.depth=[model1.depth;squeeze(-grd.z(:,ilat1,:))]; 
      model1.temp=[model1.temp;squeeze(temp1)]; 
      model1.salt=[model1.salt;squeeze(salt1)]; 
      model1.dir=romsDir{iMod}; 
      
   end
   model=[model,model1]; 
end

%% Save

save(outFile,'obs','model'); 



