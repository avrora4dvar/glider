grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow4km_r13lin_mix.nc'; 
outFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow4km_r13lin_mask.nc';

grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.mask=ncread(grdFile,'mask_rho'); 

%% plot

close all; figure; 
pcolor(grd.lon,grd.lat,grd.mask); shading flat; 
set(gca,'clim',[0,2]); 

xlim=get(gca,'xlim'); 
ylim=get(gca,'ylim'); 
d(1)=.5*diff(xlim); d(2)=.5*diff(ylim); 

%% change

while true
    [x,y,button]=ginput(1); 
    if button==1
        [dummy,iy]=min(abs(grd.lat(1,:)-y)); 
        [dummy,ix]=min(abs(grd.lon(:,iy)-x)); 
        grd.mask(ix,iy)=mod(grd.mask(ix,iy)+1,2); 
    elseif button==3
        break;
    elseif button==43
        d=.5*d; 
        xlim=x+[-1,1]*d(1); 
        ylim=y+[-1,1]*d(2);  
    elseif button==45
        d=2*d; 
        xlim=x+[-1,1]*d(1); 
        ylim=y+[-1,1]*d(2); 
    end
    
    xlim(1)=max(min(grd.lon(:)),xlim(1)); xlim(2)=min(max(grd.lon(:)),xlim(2)); 
    ylim(1)=max(min(grd.lat(:)),ylim(1)); ylim(2)=min(max(grd.lat(:)),ylim(2)); 
    pcolor(grd.lon,grd.lat,grd.mask); shading flat; 
    set(gca,'clim',[0,2],'xlim',xlim,'ylim',ylim);
end


%% save

grd.masku=grd.mask(1:end-1,:)==1 & grd.mask(2:end,:)==1; 
grd.maskv=grd.mask(:,1:end-1)==1 & grd.mask(:,2:end)==1; 

copyfile(grdFile,outFile); 
ncwrite(outFile,'mask_rho',double(grd.mask)); 
ncwrite(outFile,'mask_u',double(grd.masku)); 
ncwrite(outFile,'mask_v',double(grd.maskv)); 

