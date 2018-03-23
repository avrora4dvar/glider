clear all; clc; 
addpath('~/Matlab/figures/io'); 
addpath('~/Matlab/projects/obslist/');

winList=[2392:3:2413]; 
romsDir{1}='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_for'; 
romsDir{2}='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_for'; 
romsDir{3}='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_for'; 
romsDir{4}='/home/aruba/vol2/ipasmans/Exp/Exp35/Exp40_ana'; 

%% Read grid

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.h=ncread(grdFile,'h'); 

%% Read observations

obsDir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp37'; 
obs=[]; 
for iWin=winList
    display(iWin)
    obs1=read_obslist(fullfile(obsDir,sprintf('obslist%-.4d.nc',iWin))); 
    
    for iMod=1:length(romsDir)
        obs1.for(:,iMod)=ncread(fullfile(romsDir{iMod},sprintf('sample%-.4d.nc',iWin)),'K'); 
    end
    obs=[obs,obs1]; 
end

%% 

typeList=[4,4;6,6;7,7;8,8;1000,9999];
Mobs=zeros(size(typeList,1),1); 
M0=zeros(size(typeList,1),length(romsDir));
M1=zeros(size(M0)); 
M2=zeros(size(M0)); 
for k=1:length(obs)
    for iMod=1:length(romsDir)
        for itype=1:size(typeList,1)
            in=obs(k).type>=typeList(itype,1) & obs(k).type<=typeList(itype,2);
            in=in & obs(k).lat<48 & interp2(grd.lon',grd.lat',grd.h',obs(k).lon,obs(k).lat)<1e3; 
            Mobs(itype,1)=Mobs(itype)+sum(obs(k).sig(in).^2);
            M0(itype,iMod)=M0(itype)+sum(in);
            M1(itype,iMod)=M1(itype)+sum(obs(k).val(in)-obs(k).for(in,iMod));
            M2(itype,iMod)=M2(itype)+sum((obs(k).val(in)-obs(k).for(in,iMod)).^2);
        end
    end
end

%%

all=struct('for',[],'val',[],'type',[],'sig',[]); 
for k=1:length(obs)
    obs(k).h=interp2(grd.lon',grd.lat',grd.h',obs(k).lon(:),obs(k).lat(:)); 
    in=obs(k).lat<48 & obs(k).h<1e3; 
   
   all.sig=cat(1,all.sig,obs(k).sig(in)); 
   all.val=cat(1,all.val,obs(k).val(in)); 
   all.for=cat(1,all.for,obs(k).for(in,:)); 
   all.type=cat(1,all.type,obs(k).type(in)); 
end

%%


for iMod=1:length(romsDir)
    for itype=1:size(typeList,1)
        in=all.type>=typeList(itype,1) & all.type<=typeList(itype,2); 
        [iMod,itype,rms(all.val(in)-all.for(in,iMod)),...
            1-rms(all.val(in)-all.for(in,iMod)-mean(all.val(in)-all.for(in,iMod)))...
            /rms(abs(all.val(in)-mean(all.val(in)))+abs(all.for(in,iMod)-mean(all.for(in,iMod))))]
    end
end



