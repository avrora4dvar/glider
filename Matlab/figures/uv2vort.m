function [vort] = uv2vort(grdFile,u,v)
%uv2vort Convert velocity field to vorticity field

%% read grid

grd.lonr=ncread(grdFile,'lon_rho'); 
grd.latr=ncread(grdFile,'lat_rho'); 
grd.lonu=ncread(grdFile,'lon_u'); 
grd.latu=ncread(grdFile,'lat_u'); 
grd.lonv=ncread(grdFile,'lon_v'); 
grd.latv=ncread(grdFile,'lat_v'); 
grd.maskv=ncread(grdFile,'mask_v'); 
grd.masku=ncread(grdFile,'mask_u'); 

lon0=mean(grd.lonr(:)); lat0=mean(grd.latr(:)); 
grd.xr=earthRadius*deg2rad(grd.lonr-lon0).*cos(deg2rad(grd.latr)); 
grd.yr=earthRadius*deg2rad(grd.latr-lat0);
grd.xu=earthRadius*deg2rad(grd.lonu-lon0).*cos(deg2rad(grd.latu)); 
grd.yu=earthRadius*deg2rad(grd.latu-lat0);
grd.xv=earthRadius*deg2rad(grd.lonv-lon0).*cos(deg2rad(grd.latv)); 
grd.yv=earthRadius*deg2rad(grd.latv-lat0);

%% calculate vorticty in psi points


v(repmat(grd.maskv==0,[1,1,size(v,3)]))=0; 
u(repmat(grd.masku==0,[1,1,size(u,3)]))=0; 

dvdx=(v(2:end,:,:)-v(1:end-1,:,:))./repmat(grd.xv(2:end,:)-grd.xv(1:end-1,:),[1,1,size(v,3)]); 
dudy=(u(:,2:end,:)-u(:,1:end-1,:))./repmat(grd.yu(:,2:end)-grd.yu(:,1:end-1),[1,1,size(u,3)]); 


%% Interpolate to rho points

vort=nan(size(grd.lonr,1),size(grd.latr,2),size(u,3)); 
vort(2:end-1,2:end-1,:)=...
    .25*(dvdx(1:end-1,1:end-1,:)-dudy(1:end-1,1:end-1,:))+...
    .25*(dvdx(1:end-1,2:end,:)-dudy(1:end-1,2:end,:))+...
    .25*(dvdx(2:end,1:end-1,:)-dudy(2:end,1:end-1,:))+...
    .25*(dvdx(2:end,2:end,:)-dudy(2:end,2:end,:)); 

end

