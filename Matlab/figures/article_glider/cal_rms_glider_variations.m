%% Compare RMS different glider observations
addpath('~/Matlab/projects/obslist'); 
clear all; close all; 

romsDir={...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp38_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp39_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp41_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp42_ana'...
    }; 

obsDir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp37'; 

%% Grid

 grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
 grd.lon=ncread(grdFile,'lon_rho'); 
 grd.lat=ncread(grdFile,'lat_rho'); 
 grd.h=ncread(grdFile,'h'); 
 

%% Read observations

obs=[]; 
for iWin=2392:3:2410
   obs1=read_obslist(fullfile(obsDir,sprintf('obslist%-.d.nc',iWin)));  
   for iMod=1:length(romsDir)
       K1=ncread(fullfile(romsDir{iMod},sprintf('sample%-.4d.nc',iWin)),'K'); 
       obs1.ana(:,iMod)=K1; 
       obs1.h=interp2(grd.lon',grd.lat',grd.h',obs1.lon,obs1.lat); 
   end
   obs=[obs,obs1]; 
end

%% Combine observations

obsAll=[]; 
fnames=fields(obs1);
for iObs=1:length(obs)
    for k=1:length(fnames)
        if iObs==1; obsAll.(fnames{k})=[]; end
        obsAll.(fnames{k})=cat(1,obsAll.(fnames{k}),obs(iObs).(fnames{k})); 
    end
end

%% Calculate rms

samples.t=[2392.5:2412.5]+datenum('2005-01-01'); 
inH=obsAll.h<1e3 & obsAll.lat<48; 
for it=1:length(samples.t)
   in=obsAll.type==4 & obsAll.datenum>=samples.t(it)-.5 & obsAll.datenum<samples.t(it)+.5 & inH; 
   for k=1:size(obsAll.ana,2)
       samples.sst.mean(k,it)=mean(obsAll.ana(in,k)-obsAll.val(in)); 
       samples.sst.rms(k,it)=rms(obsAll.ana(in,k)-obsAll.val(in));
       samples.sst.corr(k,it)=corr(obsAll.ana(in,k),obsAll.val(in)); 
   end
   in=obsAll.type==6 & obsAll.datenum>=samples.t(it)-.5 & obsAll.datenum<samples.t(it)+.5 & inH; 
   for k=1:size(obsAll.ana,2)
       samples.gliderT.mean(k,it)=mean(obsAll.ana(in,k)-obsAll.val(in)); 
       samples.gliderT.rms(k,it)=rms(obsAll.ana(in,k)-obsAll.val(in));
       samples.gliderT.corr(k,it)=corr(obsAll.ana(in,k),obsAll.val(in)); 
   end
   in=obsAll.type==7 & obsAll.datenum>=samples.t(it)-.5 & obsAll.datenum<samples.t(it)+.5 & inH; 
   for k=1:size(obsAll.ana,2)
       samples.gliderS.mean(k,it)=mean(obsAll.ana(in,k)-obsAll.val(in)); 
       samples.gliderS.rms(k,it)=rms(obsAll.ana(in,k)-obsAll.val(in));
       samples.gliderS.corr(k,it)=corr(obsAll.ana(in,k),obsAll.val(in)); 
   end
   in=obsAll.type==8 & obsAll.datenum>=samples.t(it)-.5 & obsAll.datenum<samples.t(it)+.5 & inH; 
   for k=1:size(obsAll.ana,2)
       samples.uv.mean(k,it)=mean(obsAll.ana(in,k)-obsAll.val(in)); 
       samples.uv.rms(k,it)=rms(obsAll.ana(in,k)-obsAll.val(in));
       samples.uv.corr(k,it)=corr(obsAll.ana(in,k),obsAll.val(in)); 
   end
   in=obsAll.type>=1000 & obsAll.datenum>=samples.t(it)-.5 & obsAll.datenum<samples.t(it)+.5 & inH;
   if sum(in)==0
        for k=1:size(obsAll.ana,2)
           samples.ssh.mean(k,it)=NaN;
           samples.ssh.rms(k,it)=NaN; 
           samples.ssh.corr(k,it)=NaN; 
       end
   else
       for k=1:size(obsAll.ana,2)
           samples.ssh.mean(k,it)=mean(obsAll.ana(in,k)-obsAll.val(in));
           samples.ssh.rms(k,it)=rms(obsAll.ana(in,k)-obsAll.val(in));
           samples.ssh.corr(k,it)=corr(obsAll.ana(in,k),obsAll.val(in));
       end
   end


end

%% Save

save('~/Data/article_glider/glider_variation_ana.mat','samples','romsDir','obsDir'); 
