function [grdr,grdu,grdv] = read_grid_file(grdFile,nLayer)
%READ_GRID_FILE [grdr,grdu,grdv]=read_grid_file(grdFile)


grdr.lon=ncread(grdFile,'lon_rho'); 
grdr.lat=ncread(grdFile,'lat_rho'); 
grdr.mask=ncread(grdFile,'mask_rho');
grdr.h=ncread(grdFile,'h'); 
s=[-nLayer+.5:-.5]/nLayer; 
grdr.zr=roms_sigma2z(s,grdr.h); 
s=[-nLayer:0]/nLayer; 
grdr.zw=roms_sigma2z(s,grdr.h); 

grdu.lon=ncread(grdFile,'lon_u'); 
grdu.lat=ncread(grdFile,'lat_u'); 
grdu.mask=ncread(grdFile,'mask_u');
grdu.h=interp2(grdr.lon',grdr.lat',grdr.h',grdu.lon,grdu.lat); 
s=[-nLayer+.5:-.5]/nLayer; 
grdu.zr=roms_sigma2z(s,grdu.h); 
s=[-nLayer:0]/nLayer; 
grdu.zw=roms_sigma2z(s,grdu.h); 

grdv.lon=ncread(grdFile,'lon_v'); 
grdv.lat=ncread(grdFile,'lat_v'); 
grdv.mask=ncread(grdFile,'mask_v');
grdv.h=interp2(grdr.lon',grdr.lat',grdr.h',grdv.lon,grdv.lat); 
s=[-nLayer+.5:-.5]/nLayer; 
grdv.zr=roms_sigma2z(s,grdv.h); 
s=[-nLayer:0]/nLayer; 
grdv.zw=roms_sigma2z(s,grdv.h); 

end

