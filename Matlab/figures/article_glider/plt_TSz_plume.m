%% Plot z-profile of glider
clear all; clc; close all; 
addpath('~/Matlab/roms/seawater/'); 

inFile='/home/server/student/homes/ipasmans/Data/article_glider/TS_glider_expAll_for.mat';
outFile='/home/server/student/homes/ipasmans/Data/article_glider/TSz_plume.mat';
figName='/home/server/student/homes/ipasmans/Figures/article_glider/tmp'; 

lim.t=[2392,2413]+datenum('2005-01-01'); 
lim.lon=[-Inf,Inf]; 
lim.lat=[-Inf,Inf]; 
lim.z=[-Inf,Inf]; 

mode='rms'; 
field='salt'; 
zBin=[0:4:200]'; 


%% Load

load(inFile); 
model=model([4]); 

%Observations
obs.lon=repmat(obs.lon,[1,size(obs.depth,2)]); 
obs.lat=repmat(obs.lat,[1,size(obs.depth,2)]); 
obs.t=repmat(obs.t,[1,size(obs.depth,2)]); 


%% Calculate values

%Calculate average profiles 
cr=struct('t',[],'lon',[],'lat',[],'depth',[],'salt',[],'temp',[]); 
nocr=struct('t',[],'lon',[],'lat',[],'depth',[],'salt',[],'temp',[]);

%In columbia river
in=obs.salt(:,1)<=31.5 & obs.t(:,1)>=lim.t(1) & obs.t(:,1)<=lim.t(2); 
cr.t=obs.t(in,1); cr.lon=obs.lon(in,1); cr.lat=obs.lat(in,1); 
for i1=1:size(obs.salt,2)
    in1=in & ~isnan( obs.temp(:,i1) ) & ~isnan( obs.salt(:,i1) ); 
    cr.depth(i1,1)=mean(obs.depth(in1,i1)); 
    cr.salt(i1,1)=mean(obs.salt(in1,i1)); cr.salt(i1,2)=mean(model.salt(in1,i1)); 
    cr.temp(i1,1)=mean(obs.temp(in1,i1)); cr.temp(i1,2)=mean(model.temp(in1,i1)); 
end

%In columbia river
in=obs.salt(:,1)>31.5 & obs.t(:,1)>=lim.t(1) & obs.t(:,1)<=lim.t(2); 
nocr.t=obs.t(in,1); nocr.lon=obs.lon(in,1); nocr.lat=obs.lat(in,1); 
for i1=1:size(obs.salt,2)
    in1=in & ~isnan( obs.temp(:,i1) ) & ~isnan( obs.salt(:,i1) ); 
    nocr.depth(i1,1)=mean(obs.depth(in1,i1)); 
    nocr.salt(i1,1)=mean(obs.salt(in1,i1)); nocr.salt(i1,2)=mean(model.salt(in1,i1)); 
    nocr.temp(i1,1)=mean(obs.temp(in1,i1)); nocr.temp(i1,2)=mean(model.temp(in1,i1)); 
end

    
%% Save

save(outFile,'cr','nocr'); 



