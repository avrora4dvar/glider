clear all; clc; 
addpath('~/Matlab/projects/obslist'); 

obsDir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp37'; 
forDir={...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_for',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_for',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_for',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp40_ana',...
    }; 
romsDir={...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_2410/Long',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_2410/Long',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_2410/Long',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp40_ana/Long'};
romsDir=romsDir(1:end); 
%% Grid

 grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
 grd.lon=ncread(grdFile,'lon_rho'); 
 grd.lat=ncread(grdFile,'lat_rho'); 
 grd.h=ncread(grdFile,'h'); 
 
%% Forecast

samples.t=[2392:2433]; 
samples.romsDir={romsDir}; 
samples.obsDir=obsDir; 

for it=1:length(samples.t)
    if mod(samples.t(it)-samples.t(1),3)==0
        obs=read_obslist(fullfile(obsDir,sprintf('obslist%-.4d.nc',samples.t(it)))); 
        for iMod=1:length(romsDir)
            if samples.t(it)>=2416; obs.sample(:,iMod)=ncread(fullfile(romsDir{iMod},sprintf('sample%-.4d.nc',samples.t(it))),'K'); end
            if samples.t(it)<=2415; obs.sample(:,iMod)=ncread(fullfile(forDir{iMod},sprintf('sample%-.4d.nc',samples.t(it))),'K'); end
        end
        obs.h=interp2(grd.lon',grd.lat',grd.h',obs.lon,obs.lat); 
    end
    in1=obs.datenum>=samples.t(it)+datenum('2005-01-01') & obs.datenum<samples.t(it)+1+datenum('2005-01-01'); 
    in1=in1 & obs.h<1e3 & obs.lat<48; 
    for iMod=1:length(romsDir)
       in=in1 & obs.type==4; 
       samples.sst.rms(iMod,it)=rms(obs.val(in)-obs.sample(in,iMod)); 
       samples.sst.mean(iMod,it)=mean(obs.sample(in,iMod)-obs.val(in)); 
       sample.sst.corr(iMod,it)=corr( reshape(obs.val(in),[],1),reshape(obs.sample(in,iMod),[],1) ); 
       in=in1 & obs.type==6; 
       samples.gliderT.rms(iMod,it)=rms(obs.val(in)-obs.sample(in,iMod)); 
       samples.gliderT.mean(iMod,it)=mean(obs.sample(in,iMod)-obs.val(in)); 
       samples.gliderT.corr(iMod,it)=corr( reshape(obs.val(in),[],1),reshape(obs.sample(in,iMod),[],1) ); 
       in=in1 & obs.type==7; 
       samples.gliderS.rms(iMod,it)=rms(obs.val(in)-obs.sample(in,iMod)); 
       samples.gliderS.mean(iMod,it)=mean(obs.sample(in,iMod)-obs.val(in)); 
       samples.gliderS.corr(iMod,it)=corr( reshape(obs.val(in),[],1),reshape(obs.sample(in,iMod),[],1) ); 
       in=in1 & obs.type==8; 
       samples.uv.rms(iMod,it)=rms(obs.val(in)-obs.sample(in,iMod)); 
       samples.uv.mean(iMod,it)=mean(obs.sample(in,iMod)-obs.val(in)); 
       samples.uv.corr(iMod,it)=corr( reshape(obs.val(in),[],1),reshape(obs.sample(in,iMod),[],1) ); 
       in=in1 & obs.type>=1e3; 
       if sum(in)==0
           samples.ssh.rms(iMod,it)=NaN;
           samples.ssh.mean(iMod,it)=NaN;
           samples.ssh.corr(iMod,it)=NaN;
       else
           samples.ssh.rms(iMod,it)=rms(obs.val(in)-obs.sample(in,iMod));
           samples.ssh.mean(iMod,it)=mean(obs.sample(in,iMod)-obs.val(in));
           samples.ssh.corr(iMod,it)=corr( reshape(obs.val(in),[],1),reshape(obs.sample(in,iMod),[],1) );
       end
    end
end

samples.t=samples.t+datenum('2005-01-01'); 
save('~/Data/article_glider/for_long_sample.mat','samples','grd'); 
%% 

cmap=lines(length(romsDir)); 
if false
    close all;
    fnames={'sst','gliderT','gliderS','uv','ssh'};
    for iF=1:length(fnames)
        figure(); hold on;
        for iMod=1:length(romsDir)
            plot(samples.t,samples.(fnames{iF}).rms(iMod,:),'color',cmap(iMod,:));
        end
        set(gca,'xlim',[2410,2434]+datenum('2005-01-01')); 
        datetick('x'); 
        title(fnames{iF});
        legend({'exp35','exp36','exp37','free'});
        
    end
    
    
end
    
