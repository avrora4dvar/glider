function [ output_args ] = line2index( input_args )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

%% read grid

grdr.lon=ncread(grid_file,'lon_rho'); 
grdr.lat=ncread(grid_file,'lat_rho'); 
grdr.mask=ncread(grid_file,'mask_rho');
grdr.h=ncread(grid_file,'h'); 

deg2m=[cos(deg2rad(mean(grdr.lat(:)))),1]*2*pi*earthRadius;
degcen=[mean(grdr.lon(:)),mean(grdr.lat(:))]; 

grdp.lon=.25*( grdr.lon(1:end-1,1:end-1)+grdr.lon(2:end,1:end-1)+...
    grdr.lon(1:end-1,2:end)+grdr.lon(2:end,2:end)); 
grdp.lat=.25*( grdr.lat(1:end-1,1:end-1)+grdr.lat(2:end,1:end-1)+...
    grdr.lat(1:end-1,2:end)+grdr.lat(2:end,2:end)); 
grdp.dlon=grdp.lon(2:end,:)-grdp.lon(1:end-1,:); 
grdp.dlat=grdp.lat(:,2:end)-grdp.lat(:,1:end-1); 
grdp.x=(grdp.lon-degcen(1))*deg2m(1); 
grdp.y=(grdp.lat-degcen(2))*deg2m(2); 

grdu.lon=ncread(grid_file,'lon_u'); 
grdu.lat=ncread(grid_file,'lat_u'); 
grdu.mask=ncread(grid_file,'mask_u');
grdu.h=.5*grdr.h(1:end-1,:)+.5*grdr.h(2:end,:); 

grdv.lon=ncread(grid_file,'lon_v'); 
grdv.lat=ncread(grid_file,'lat_v'); 
grdv.mask=ncread(grid_file,'mask_v');
grdv.h=.5*grdr.h(:,1:end-1)+.5*grdr.h(:,2:end); 

dl=min(min(grdp.dlon(:)*deg2m(1)),min(grdp.dlat(:)*deg2m(2))); 

%% find u,v points along the line

sec.x=deg2m(1)*(loc.lon-degcen(1)); 
sec.y=deg2m(2)*(loc.lat-degcen(2)); 
sec.d=[0,hypot(diff(sec.x),diff(sec.y))]; 

sec.di=[0:.25*dl:max(sec.d)-1e-6,max(sec.d)]; 
sec.xi=interp1(sec.d,sec.x,sec.di); 
sec.yi=interp1(sec.d,sec.y,sec.di); 

for i0=1:length(sec.di)
   [dummy,sec.ilonlat(i0,1)]=min( abs(grdp.x(:,1)-sec.xi(i0)) ); 
   [dummy,sec.ilonlat(i0,2)]=min( abs(grdp.y(1,:)-sec.yi(i0)) ); 
end

[uniLL,LL2uni,uni2LL]=unique(sec.ilonlat,'rows','first');

if size(uniLL,1)<=1
    error('Insufficient points to make a section'); 
else
    sec.ilon=[]; sec.ilat=[]; sec.dim=[]; 
    for i0=2:size(uniLL,1)
        stepdiff=uniLL(i0,:)-uniLL(i0-1,:);       
       if stepdiff(1)==-1 
           sec.dim=[sec.dim,2]; 
           sec.ilon=[sec.ilon,uniLL(i0,1)]; 
           sec.ilat=[sec.ilat,uniLL(i0,2)];
       end
       if stepdiff(1)==1 && uniLL(i0,1)+1<=size(grdv.lon,1)
           sec.dim=[sec.dim,-2]; 
           sec.ilon=[sec.ilon,uniLL(i0,1)+1]; 
           sec.ilat=[sec.ilat,uniLL(i0,2)]; 
       end
       uniLL(i0,1)=uniLL(i0,1)+stepdiff(1); 
       if stepdiff(2)==-1
           sec.dim=[sec.dim,-1]; 
           sec.ilon=[sec.ilon,uniLL(i0,1)]; 
           sec.ilat=[sec.ilat,uniLL(i0,2)]; 
       end
       if stepdiff(2)==1 && uniLL(i0,2)+1<=size(grdu.lat,2)
           sec.dim=[sec.dim,1];
           sec.ilon=[sec.ilon,uniLL(i0,1)]; 
           sec.ilat=[sec.ilat,uniLL(i0,2)+1]; 
       end
       if any(abs(stepdiff)>1)
           error('point has been skipped. Use smaller dl.'); 
       end
    end
    
end

end

