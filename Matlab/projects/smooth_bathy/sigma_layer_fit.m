%% Bathymetry interpolation that places sigma-points at same z
clear all; clc; 
addpath('~/Matlab/figures/'); 


grdFileIn='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc';
grdFileOut='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow4km_r13lin_mix.nc';
outFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow4km_r13lin_fit.nc';
sigma=[-39.5:-.5]/40; 

%% read grid

grd1.lon=ncread(grdFileIn,'lon_rho'); 
grd1.lat=ncread(grdFileIn,'lat_rho'); 
grd1.mask=ncread(grdFileIn,'mask_rho'); 
grd1.h=ncread(grdFileIn,'h'); 

grd2.lon=ncread(grdFileOut,'lon_rho'); 
grd2.lat=ncread(grdFileOut,'lat_rho'); 
grd2.mask=ncread(grdFileOut,'mask_rho'); 

%% Initial interpolation

grd1.z=roms_sigma2z(sigma,grd1.h,0*grd1.h); 
for i3=1:size(grd1.z,3)
    grd2.zi(:,:,i3)=interp2(grd1.lon',grd1.lat',squeeze(grd1.z(:,:,i3))',grd2.lon,grd2.lat); 
end
grd2.hi=interp2(grd1.lon',grd1.lat',squeeze(grd1.h)',grd2.lon,grd2.lat); 
grd2.maski=interp2(grd1.lon',grd1.lat',squeeze(grd1.mask)',grd2.lon,grd2.lat); 

%% Start fitting

grd2.h=nan(size(grd2.lon)); 
grd2.J=nan(size(grd2.lon)); 
for i1=1:size(grd2.lon,1)
    display(i1); 
    for i2=1:size(grd2.lat,2)
        
        if grd2.mask(i1,i2)==0; continue; end
        z0=squeeze(grd2.zi(i1,i2,:)); 
        
        h=[grd2.hi(i1,i2)]; 
        z1=squeeze(roms_sigma2z(sigma,h(end),0)); 
        J=[max(abs(z1-z0))]; 
        
        %Find depth for difference is minimal
        dir=1; step=.05*h(1); 
        for k=1:500
            h=[h,h(end)+step*dir]; 
            z1=squeeze(roms_sigma2z(sigma,h(end),0)); 
            J=[J,max(abs(z1-z0))]; 
            if sign(J(end)-J(end-1))>0; dir=-dir; step=.5*step; end
            if abs(J(end)-J(end-1))<.05*J(1) && J(end)<J(1); break; end
        end
        [minval,minloc]=min(J); 
        grd2.h(i1,i2)=h(minloc); 
        grd2.J(i1,i2)=J(minloc);   
       

    end
end

%Apply mask
grd2.h(grd2.mask==0)=3; 

%% Save to file

copyfile(grdFileOut,outFile);
grd2.zr=roms_sigma2z([-39.5:-.5]/40,grd2.h,0*grd2.h); 
grd2.zw=roms_sigma2z([-40:0]/40,grd2.h,0*grd2.h); 
ncwrite(outFile,'h',grd2.h); 
ncwrite(outFile,'z0_r',grd2.zr); 
ncwrite(outFile,'z0_w',grd2.zw); 
