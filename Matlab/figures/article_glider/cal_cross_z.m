%calculate a zonal cross section
clear all; clc; close all; 

addpath('/home/server/student/homes/ipasmans/Matlab/roms/seawater/'); 
addpath('..'); 

%% input

%grid
grid_file='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
%roms dir
roms_dir={'/home/aruba/vol2/ipasmans/Exp/Exp35/Exp35_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp36_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp37_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp38_ana',...
    '/home/aruba/vol2/ipasmans/Exp/Exp35/Exp40_ana',...
    }; 
%output file
out_file='/home/server/student/homes/ipasmans/Data/article_glider/cross_NH10_his_zrho_ana.mat'; 
%latitude cross-section
loc.lat=43.7; %NH10 line=44.6 - 44.7->44.64
loc.lat=44.64;
loc.lon=[-129,-123.5]; 
loc.t=[datenum('2011-04-01 00:00'):1:datenum('2011-05-27 00:00')]; 
loc.t=[datenum('2011-07-21 00:00'):1/24:datenum('2011-08-14 00:00')]; 
loc.t=[2392:1/24:2413]+datenum('2005-01-01'); 

%% read grid

%rho-grid
grd.rho.lon=ncread(grid_file,'lon_rho'); 
grd.rho.lat=ncread(grid_file,'lat_rho'); 
grd.rho.mask=ncread(grid_file,'mask_rho'); 
grd.rho.z=ncread(grid_file,'z0_r'); 
grd.rho.dz=diff(ncread(grid_file,'z0_w'),[],3); 

%u-grid
grd.u.lon=ncread(grid_file,'lon_u'); 
grd.u.lat=ncread(grid_file,'lat_u'); 
grd.u.mask=ncread(grid_file,'mask_u'); 

%v-grid
grd.v.lon=ncread(grid_file,'lon_v'); 
grd.v.lat=ncread(grid_file,'lat_v'); 
grd.v.mask=ncread(grid_file,'mask_v'); 

%% Calculate grid

w=interp1(grd.rho.lat(1,:),[1:size(grd.rho.lat,2)],loc.lat); 
in=[find(grd.rho.lon(:,1)<min(loc.lon),1,'last'):find(grd.rho.lon(:,1)>max(loc.lon),1,'first')]; 
z=(w-floor(w))*grd.rho.z(in,ceil(w),:)+(ceil(w)-w)*grd.rho.z(in,floor(w),:); z=squeeze(z); 
dz=(w-floor(w))*grd.rho.dz(in,ceil(w),:)+(ceil(w)-w)*grd.rho.dz(in,floor(w),:); dz=squeeze(dz);
p1=sw_pres(-z,loc.lat); 

%% read model

levels=[1020:.5:1030]; 
for iMod=1:length(roms_dir)
   display(roms_dir{iMod}); 
    
   %temperature
   temp=get_roms_cross_zonal(grid_file,roms_dir{iMod},'temp',grd.rho,loc); 
   %salinity
   salt=get_roms_cross_zonal(grid_file,roms_dir{iMod},'salt',grd.rho,loc);
   %rho
   for it=1:length(temp.t)
       display(datestr(temp.t(it),'yyyymmdd HH:MM')); 
       rho1=sw_pden(squeeze(salt.val(:,2:end-1,it)),squeeze(temp.val(:,2:end-1,it)),p1,0); 
       rho1=stable_rho(reshape(rho1,[size(rho1,1),1,size(rho1,2)]),...
           reshape(dz,[size(rho1,1),1,size(rho1,2)]),'min_step',1e-6);
       rho1=squeeze(rho1); 
       for i1=1:size(rho1,1)
           if sum(~isnan(rho1(i1,:)))>2 
            z1(i1,:,it)=interp1(rho1(i1,:),z(i1,:),levels); 
           else
               z1(i1,:,it)=NaN; 
           end
       end
   end

   
   
   %save output
   model(iMod)=struct('t',temp.t,'lon',temp.lon,'z',z1,...
       'roms_dir',roms_dir{iMod},'rho',levels); 
end

save(out_file,'model','loc'); 




