%% Get cross-section with forecast/analysis
addpath('/home/server/student/homes/ipasmans/Matlab/roms/seawater');
addpath('/home/server/student/homes/ipasmans/Matlab/projects/obslist/'); 
clear all; clc; 

obsFile='/home/aruba/vol2/ipasmans/Exp/Obs/Exp36/obslist2392.nc'; 
grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
expDir='/home/aruba/vol2/ipasmans/Exp/Exp35'; 
expName={'Exp35','Exp36','Exp37'}; 
outFile='/home/server/student/homes/ipasmans/Data/article_glider/cross_anafor.mat'; 
dateRef=datenum('2005-01-01'); 


%% Read grid

grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.mask=ncread(grdFile,'mask_rho'); 
grd.zr=ncread(grdFile,'z0_r'); 
grd.zw=ncread(grdFile,'z0_w'); 

p=sw_pres(-grd.zr(:),grd.lat(1)); 

%% Observations

obs=read_obslist(obsFile); 
[~,ilat]=min(abs(mean(obs.lat(obs.type==6))-grd.lat(1,:))); 
error('ok')

grd.lon=grd.lon(:,ilat); 
grd.lat=grd.lat(:,ilat); 
grd.mask=grd.mask(:,ilat); 
grd.zr=squeeze(grd.zr(:,ilat,:)); 
grd.zw=squeeze(grd.zw(:,ilat,:)); 


%% Read file

for iExp=1:length(expName)
    model(iExp)=struct(...
        't',[],...
        'temp',struct('for',[],'ana',[]),...
        'salt',struct('for',[],'ana',[]),...
        'rho',struct('for',[],'ana',[]),...
        'v',struct('for',[],'ana',[]) ); 
end

winNo=[2392:3:2410]; 
for iWin=winNo
   for iExp=1:length(expName)
      %Forecast
      fname=fullfile(expDir,sprintf('%s_for',expName{iExp}),sprintf('ocean_avg_%-.4d.nc',iWin)); 
      display(fname); 
      t1=ncread(fname,'ocean_time')/3600/24+dateRef; 
      model(iExp).t=cat(1,model(iExp).t,t1); 
      valTemp=ncread(fname,'temp',[1,ilat,1,1],[Inf,1,Inf,Inf]); valTemp=squeeze(valTemp); 
      model(iExp).temp.for=cat(3,model(iExp).temp.for,valTemp); 
      valSalt=ncread(fname,'salt',[1,ilat,1,1],[Inf,1,Inf,Inf]); valSalt=squeeze(valSalt); 
      model(iExp).salt.for=cat(3,model(iExp).salt.for,valSalt); 
      p1=repmat(p(:),[size(valTemp,3),1]); 
      rho=sw_pden(valSalt(:),valTemp(:),p1,0); rho=reshape(rho,size(valTemp)); 
      model(iExp).rho.for=cat(3,model(iExp).rho.for,rho); 
      valV=ncread(fname,'v',[1,ilat-1,1,1],[Inf,2,Inf,Inf]); valV=squeeze(mean(valV,2)); 
      model(iExp).v.for=cat(3,model(iExp).v.for,valV); 
      
      %Analysis
      fname=fullfile(expDir,sprintf('%s_ana',expName{iExp}),sprintf('ocean_avg_%-.4d.nc',iWin));
      display(fname); 
      valTemp=ncread(fname,'temp',[1,ilat,1,1],[Inf,1,Inf,Inf]); valTemp=squeeze(valTemp); 
      model(iExp).temp.ana=cat(3,model(iExp).temp.ana,valTemp); 
      valSalt=ncread(fname,'salt',[1,ilat,1,1],[Inf,1,Inf,Inf]); valSalt=squeeze(valSalt); 
      model(iExp).salt.ana=cat(3,model(iExp).salt.ana,valSalt); 
      p1=repmat(p(:),[size(valTemp,3),1]); 
      rho=sw_pden(valSalt(:),valTemp(:),p1,0); rho=reshape(rho,size(valTemp)); 
      model(iExp).rho.ana=cat(3,model(iExp).rho.ana,rho); 
      valV=ncread(fname,'v',[1,ilat-1,1,1],[Inf,2,Inf,Inf]); valV=squeeze(mean(valV,2)); 
      model(iExp).v.ana=cat(3,model(iExp).v.ana,valV); 
      
      
   end
end

%% Save

save(outFile,'model','expName','expDir','grd','ilat'); 

