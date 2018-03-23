%generate bathymetry
clear all; clc; close all; 

outFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow4km_r18lin.nc'; 
%% load grid

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow4km_roms.nc';
grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.h=ncread(grdFile,'h'); 
grd.mask=ncread(grdFile,'mask_rho')==1; 
grd.h(~grd.mask)=NaN; 
hc=3;

%% load h

load('bathymetry_com4km.mat'); 
h=max(h,10); 
%h(~grd.mask)=3; 

%% Add Columbia river

% h(237:275,308:318)=15; 
% h(~grd.mask)=hc; 

h0=h; 
h0(~grd.mask)=NaN; 

%plot 
figure(); 
pcolor(grd.lon,grd.lat,-grd.h); shading flat; set(gca,'clim',[-4e3,0]); 
colorbar; 

%plot 
figure(); 
pcolor(grd.lon,grd.lat,-h0); shading flat; set(gca,'clim',[-4e3,0]); 
colorbar; 

%% Calculate rx0/rx1

%rx0
rx0=beckman_haidvogel_haney(h0); 

%rx1
sigma=[-39.5:-.5]/40;
z=roms_sigma2z(sigma,h0); 
rx1=beckman_haidvogel_haney(z); 

display(sprintf('min/max rx0: %.2f %.2f',[min(rx0(:)),max(rx0(:))])); 
display(sprintf('min/max rx1: %.2f %.2f',[min(rx1(:)),max(rx1(:))])); 

%plot
figure(); 
pcolor(grd.lon,grd.lat,rx1); shading flat; 
colorbar; 

%% smooth bathymetry

hs=smooth_meo(h0,grd.mask,.18);
hs(isnan(hs))=10; 
hs=smooth_lin(h,hs,grd.mask);
hs(~grd.mask)=NaN; 

%plot 
figure(); 
pcolor(grd.lon,grd.lat,-hs); shading flat; set(gca,'clim',[-4e3,0]); 
colorbar; 

%plot 
figure(); 
pcolor(grd.lon,grd.lat,-hs+h0); shading flat; set(gca,'clim',[-500,500]); 
colorbar; 
%% Calculate rx0/rx1

%rx0
rx0=beckman_haidvogel_haney(hs); 

%rx1
sigma=[-39.5:-.5]/40;
z=roms_sigma2z(sigma,hs); 
rx1=beckman_haidvogel_haney(z); 

display(sprintf('min/max rx0: %.2f %.2f',[min(rx0(:)),max(rx0(:))])); 
display(sprintf('min/max rx1: %.2f %.2f',[min(rx1(:)),max(rx1(:))])); 


%% write

%nan-value
hs(isnan(hs))=hc; 
hs=max(hs,hc); 

copyfile(grdFile,outFile); 
ncwrite(outFile,'h',hs); 

sigma=[-39.5:-.5]/40;
z=roms_sigma2z(sigma,hs); 
ncwrite(outFile,'z0_r',z); 

sigma=[-40:0]/40;
z=roms_sigma2z(sigma,hs); 
ncwrite(outFile,'z0_w',z); 


