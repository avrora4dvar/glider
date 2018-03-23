%% make new bathymetry for Salish Sea
clear all; clc; 

out_file='bathymetry_com4km.mat'; 

%% read grid

grdFile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow4km_r13lin_mix.nc';

display(sprintf('read %s',grdFile)); 
grd.lon=ncread(grdFile,'lon_rho'); 
grd.lat=ncread(grdFile,'lat_rho'); 
grd.mask=ncread(grdFile,'mask_rho')==1; 
grd.h=ncread(grdFile,'h'); 
deg2m=2*pi*earthRadius/360*[cos(deg2rad(mean(grd.lat(:)))),1]; 
 
%% read ETOPO

topFile='/home/muikku/kurapov/BATHYMETRY_OR_12sec/uswest_etopo1.nc';
display(sprintf('read %s',topFile)); 
topo.lon=ncread(topFile,'lon'); 
topo.lat=ncread(topFile,'lat'); 
topo.h=-ncread(topFile,'z'); 

%interpolate
topo.hi=interp2(topo.lon',topo.lat',topo.h',grd.lon,grd.lat); 
topo.maski=interp2(topo.lon',topo.lat',double(topo.h'>3),grd.lon,grd.lat); 
topo.maski=topo.maski==1; 

%% read 12 sec bathmetry

display('read bathy12s'); 
bathy12s=load('/home/muikku/kurapov/BATHYMETRY_OR_12sec/bathy_12secXYZ.mat'); 
bathy12s.xbathy=repmat(bathy12s.xbathy,[size(bathy12s.zbathy,1),1]); 
bathy12s.ybathy=repmat(bathy12s.ybathy,[1,size(bathy12s.zbathy,2)]); 

%
bnd1=.5*grd.lon(1:end-1,1)+.5*grd.lon(2:end,1); 
bnd2=.5*grd.lat(1,1:end-1)+.5*grd.lat(1,2:end); 
bathy12s.hi=nan(size(grd.lon,1),size(grd.lon,2)); 
bathy12s.maski=zeros(size(grd.lon,1),size(grd.lon,2)); 
for i1=1:length(bnd1)-1
    in1=bathy12s.xbathy(1,:)>=bnd1(i1) & bathy12s.xbathy(1,:)<=bnd1(i1+1); 
    if ~any(in1); continue; end
    for i2=1:length(bnd2)-1
        in2=bathy12s.ybathy(:,1)>=bnd2(i2) & bathy12s.ybathy(:,1)<=bnd2(i2+1); 
        if ~any(in2); continue; end
        
        tmp=bathy12s.zbathy(in2,in1); 
        if sum(tmp(:)<3) >.5*length(tmp(:))
            bathy12s.hi(i1+1,i2+1)=3; 
            bathy12s.maski(i1+1,i2+1)=0; 
        else
            tmp(~isnan(tmp) & tmp<3)=3; 
            bathy12s.hi(i1+1,i2+1)=nanmean(tmp(:)); 
            bathy12s.maski(i1+1,i2+1)=sum(~isnan(tmp(:)));  
        end
    end    
end

error('ok'); 
%% associate points ETOPO

display('associate ETOPO'); 

if true
%for each ETOPO point find closest grid point
topo.ilon=nan(size(topo.lon)); 
topo.ilat=nan(size(topo.lat)); 

%remove points outside grid
inan= topo.lon<=max(grd.lon(:)) & topo.lon>=min(grd.lon(:));
inan= inan & topo.lat<=max(grd.lat(:)) & topo.lat>=min(grd.lat(:));

%output grid
topo.hmax=-Inf*ones(size(grd.mask)); 
topo.hmean=zeros(size(grd.mask)); 
topo.hn=zeros(size(grd.mask)); 

for i2=1:size(topo.lon,2)
    for i1=1:size(topo.lat,1)
        if ~inan(i1,i2); continue; end
        d=(grd.lon(:,1)-topo.lon(i1,i2)); 
        [dummy,j1]=min(abs(d)); 
        d=(grd.lat(1,:)-topo.lat(i1,i2)); 
        [dummy,j2]=min(abs(d)); 
%         d=hypot(deg2m(1)*(grd.lon-topo.lon(i1,i2)),deg2m(2)*(grd.lat-topo.lat(i1,i2))); 
%         [j1,j2]=find(d==min(d(:)),1); 

        %output
        topo.hn(j1,j2)=topo.hn(j1,j2)+1;
        topo.hmax(j1,j2)=max(topo.hmax(j1,j2),topo.h(i1,i2)); 
        topo.hmean(j1,j2)=topo.hmean(j1,j2)+topo.h(i1,i2); 
    end
end

topo.hmax(topo.hn==0)=NaN; 
topo.hmean(topo.hn==0)=NaN; 
topo.hmean=topo.hmean./topo.hn; 

else
   topo.hinterp=interp2(topo.lon',topo.lat',topo.h',grd.lon,grd.lat);  
end

%% associate points BATHY12s

display('associate BATHY12s'); 

%for each ETOPO point find closest grid point
bathy12s.ilon=nan(size(bathy12s.ybathy)); 
bathy12s.ilat=nan(size(bathy12s.xbathy)); 

%remove points outside grid
inan= bathy12s.xbathy<=max(grd.lon(:)) & bathy12s.xbathy>=min(grd.lon(:));
inan= inan & bathy12s.ybathy<=max(grd.lat(:)) & bathy12s.ybathy>=min(grd.lat(:));

%output grid
bathy12s.hn=zeros(size(grd.mask)); 
bathy12s.hmean=zeros(size(grd.mask)); 
bathy12s.hmax=-Inf*ones(size(grd.mask)); 
bathy12s.hmin=Inf*ones(size(grd.mask)); 

for i2=1:size(bathy12s.xbathy,2)
    for i1=1:size(bathy12s.xbathy,1)
        if ~inan(i1,i2); continue; end
        d=grd.lon(:,1)-bathy12s.xbathy(i1,i2); 
        [dummy,j1]=min(abs(d)); 
        d=grd.lat(1,:)-bathy12s.ybathy(i1,i2); 
        [dummy,j2]=min(abs(d)); 
%         d=hypot(deg2m(1)*(grd.lon-bathy12s.xbathy(i1,i2)),deg2m(2)*(grd.lat-bathy12s.ybathy(i1,i2))); 
%         [j1,j2]=find(d==min(d(:)),1); 
        bathy12s.ilon(i1,i2)=j1; bathy12s.ilat(i1,i2)=j2; 
        
        
        %output
        bathy12s.hn(j1,j2)=bathy12s.hn(j1,j2)+1;
        bathy12s.hmean(j1,j2)=bathy12s.hmean(j1,j2)+bathy12s.zbathy(i1,i2);
        bathy12s.hmin(j1,j2)=min(bathy12s.hmin(j1,j2),bathy12s.zbathy(i1,i2)); 
        bathy12s.hmax(j1,j2)=max(bathy12s.hmax(j1,j2),bathy12s.zbathy(i1,i2)); 
    end
end

bathy12s.hmean(bathy12s.hn==0)=NaN; 
bathy12s.hmax(bathy12s.hn==0)=NaN; 
bathy12s.hmin(bathy12s.hn==0)=NaN; 
bathy12s.hmean=bathy12s.hmean./bathy12s.hn; 
bathy12s.hinterp=interp2(bathy12s.xbathy,bathy12s.ybathy,bathy12s.zbathy,grd.lon,grd.lat); 
bathy12s.hinterp(bathy12s.hn==0)=NaN; 



%% select whether to use maximal bathymetry or mean bathymetry

display('select what bathymetry to use'); 

h=nan(size(grd.mask)); hmean=h; hmax=h; 

inan=bathy12s.hmean<100; 
h(inan)=bathy12s.hmin(inan); 
inan=bathy12s.hmean>=100;
h(inan)=bathy12s.hmean(inan); 
inan=isnan(bathy12s.hmean); 
h(inan)=topo.hmean(inan); 

sum(isnan(h))
%% interpolate/extrapolate

display('interpolating/extrapolating'); 
while any(isnan(h))
   
  for i2=1:size(h,2)
      for i1=1:size(h,1)
          if isnan(h(i1,i2))
              hloc= h(max(1,i1-1):min(size(h,1),i1+1),...
                  max(1,i2-1):min(size(h,2),i2+1));
              h(i1,i2)=nanmean(hloc(:)); 
          end
      end
  end
              
end
        
%% save

save(out_file,'h','grdFile'); 

