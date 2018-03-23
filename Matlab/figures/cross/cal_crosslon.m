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
out_file='/home/server/student/homes/ipasmans/Data/article_glider/cross_NH10_his_ana.mat'; 
%latitude cross-section
loc.lat=43.7; %NH10 line=44.6 - 44.7->44.64
loc.lat=44.64;
loc.lon=[-129,-123.5]; 
loc.t=[datenum('2011-04-01 00:00'):1:datenum('2011-05-27 00:00')]; 
loc.t=[datenum('2011-07-21 00:00'):2/24:datenum('2011-08-14 00:00')]; 

%% read grid

%rho-grid
grd.rho.lon=ncread(grid_file,'lon_rho'); 
grd.rho.lat=ncread(grid_file,'lat_rho'); 
grd.rho.mask=ncread(grid_file,'mask_rho'); 

%u-grid
grd.u.lon=ncread(grid_file,'lon_u'); 
grd.u.lat=ncread(grid_file,'lat_u'); 
grd.u.mask=ncread(grid_file,'mask_u'); 

%v-grid
grd.v.lon=ncread(grid_file,'lon_v'); 
grd.v.lat=ncread(grid_file,'lat_v'); 
grd.v.mask=ncread(grid_file,'mask_v'); 

%% read model

for iMod=1:length(roms_dir)
   display(roms_dir{iMod}); 
    
   %temperature
   temp=get_roms_cross_zonal(grid_file,roms_dir{iMod},'temp',grd.rho,loc); 
   %salinity
   salt=get_roms_cross_zonal(grid_file,roms_dir{iMod},'salt',grd.rho,loc);
   %v
   v=get_roms_cross_zonal(grid_file,roms_dir{iMod},'v',grd.v,loc);
   %dudx
   dudy=get_roms_cross_zonal(grid_file,roms_dir{iMod},'u_dy',grd.u,loc);
   
   %rel vorticity
   dx=deg2rad(v.lon(2:end,1)-v.lon(1:end-1,1))*earthRadius*cos(deg2rad(loc.lat)); 
   dvdx=(v.val(2:end,:,:)-v.val(1:end-1,:,:))./repmat(dx(:),[1,size(v.val,2),size(v.val,3)]); 
   dvdx=interp1(.5*v.lon(2:end,1)+.5*v.lon(1:end-1,1),dvdx,dudy.lon(:,1)); 
   vort=dvdx-dudy.val; vort=interp1(dudy.lon(:,1),vort,temp.lon(:,1)); 
   
   %save output
   model(iMod)=struct('t',temp.t,'lon',temp.lon,'z',temp.z,...
       'temp',temp.val,'salt',salt.val,...
       'vort',vort,'roms_dir',roms_dir{iMod}); 
end

save(out_file,'model','loc'); 




