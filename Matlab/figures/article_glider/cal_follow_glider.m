clear all; clc; 
addpath('/home/server/student/homes/ipasmans/Matlab/figures/io'); 
addpath('/home/server/student/homes/ipasmans/Matlab/figures/'); 

romsDir={...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp40_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_ana',...
    };  

%% Read grid

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.zr=ncread(grdFile,'z0_r'); 
grd.zw=ncread(grdFile,'z0_w'); 
grd.h=ncread(grdFile,'h'); 

z=[2:4:198]; 
%% Read obs

obsDir='/home/server/student/homes/ipasmans/Obs/glider/osu'; 
obs1=io_glider_osu(obsDir,'dir'); 

in=obs1.t>=datenum('2011-07-21') & obs1.t<=datenum('2011-08-11'); 
obs1.lon=obs1.lon(in); obs1.lat=obs1.lat(in); obs1.depth=obs1.depth(in,:); obs1.t=obs1.t(in); 
obs1.temp=obs1.temp(in,:); obs1.salt=obs1.salt(in,:); 


obs.depth=reshape([2:4:198],1,[]); 
obs.t=obs1.t; obs.lon=obs1.lon; obs.lat=obs1.lat; 
obs.temp=nan(length(obs1.t),length(obs.depth)); 
obs.salt=nan(length(obs1.t),length(obs.depth));
for i1=1:size(obs1.temp,1)
    in=~isnan(obs1.depth(i1,:)); 
    obs.temp(i1,:)=interp1(obs1.depth(i1,in),obs1.temp(i1,in),obs.depth); 
    obs.salt(i1,:)=interp1(obs1.depth(i1,in),obs1.salt(i1,in),obs.depth); 
end

%% Read model

model=[]; 
for iMod=1:length(romsDir)
    
    model1.temp=nan(size(obs.temp)); 
    model1.salt=nan(size(obs.salt)); 
        
    loc.t=repmat(obs.t,[1,size(obs.temp,2)]); 
    loc.lon=repmat(obs.lon,[1,size(obs.temp,2)]); 
    loc.lat=repmat(obs.lat,[1,size(obs.temp,2)]); 
    loc.depth=repmat(obs.depth,[size(obs.temp,1),1]); 
    
    loc.t=reshape(loc.t,1,[]); 
    loc.lon=reshape(loc.lon,1,[]); 
    loc.lat=reshape(loc.lat,1,[]); 
    loc.depth=reshape(loc.depth,1,[]); 
    
    display(sprintf('%s-temp',romsDir{iMod})); 
    out1=get_roms_3d(grdFile,romsDir{iMod},'temp',grd,loc); 
    display(sprintf('%s-salt',romsDir{iMod})); 
    out2=get_roms_3d(grdFile,romsDir{iMod},'salt',grd,loc); 
    
    model1.temp=reshape(out1.val,size(obs.temp)); 
    model1.salt=reshape(out2.val,size(obs.salt)); 
    model1.z=reshape(out1.z,size(obs.salt)); 
    model1.dir=romsDir{iMod}; 
    
    %Verwijder extrapolatie
    h1=interp2(grd.lon',grd.lat',grd.h',obs.lon,obs.lat); 
    for i1=1:length(h1)
        in=obs.depth>=h1(i1); 
        model1.temp(i1,in)=NaN; 
        model1.salt(i1,in)=NaN; 
    end
    
    model=[model,model1]; 
    
    
end

%% Save

save('/home/server/student/homes/ipasmans/Data/article_glider/along_glider.mat','obs','model'); 
