grid_file='/home/aruba/vol1/ipasmans/stampede/Prm/grd_ow2km_roms.nc';
out_file='/home/server/student/homes/ipasmans/PhD/projects/boundary/bnd_test.nc'; 
tList=datenum({'2011-04-01'}); 

%% read output grid file

%x:longitudes 0-360
%y:latitudes

%read grid coordinates output grid
grd.x=ncread(grid_file,'lon_rho'); grd.x=mod(grd.x,360); 
grd.y=ncread(grid_file,'lat_rho'); 
grd.mask=ncread(grid_file,'mask_rho'); 
gsize=[size(grd.x),length(tList)]; 

%read grid coordinate input grid
grib.grid=dlmread('/home/server/student/homes/ipasmans/PhD/projects/boundary/NAM/latlon-g218.txt'); 
bsize=[max(grib.grid(:,1)),max(grib.grid(:,2))]; 
grib.x=mod(-reshape(grib.grid(:,4),bsize),360); 
grib.y=reshape(grib.grid(:,3),bsize); 

%select subset of input grid
grib.mask1=[Inf,-Inf]; grib.mask2=[Inf,-Inf]; 
for i2=1:size(grib.x,2)
    grib.mask1(1)=min(grib.mask1(1),find(grib.x(:,i2)>=min(grd.x(:)),1,'first')-1); 
    grib.mask1(2)=max(grib.mask1(2),find(grib.x(:,i2)<=max(grd.x(:)),1,'last')+1); 
end
for i1=1:size(grib.x,1)
    grib.mask2(1)=min(grib.mask2(1),find(grib.y(i1,:)>=min(grd.y(:)),1,'first')-1); 
    grib.mask2(2)=max(grib.mask2(2),find(grib.y(i1,:)<=max(grd.y(:)),1,'last')+1); 
end
grib.mask1=[grib.mask1(1):grib.mask1(2)]; 
grib.mask2=[grib.mask2(1):grib.mask2(2)]; 

%% create output file

%create output file
if exist(out_file,'file')
    delete(out_file); 
end
varname='Tair'; 
nccreate(out_file,sprintf('%s_time',varname),'Format','netcdf4_classic',...
    'Dimensions',{'time',Inf}); 
nccreate(out_file,sprintf('%s',varname),'Format','netcdf4_classic',...
    'Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf}); 

%write time
ncwrite(out_file,sprintf('%s_time',varname),tList-datenum('2005-01-01')); 

%% 

addpath('/home/server/student/homes/ipasmans/OpenEarth/nctoolbox'); 
setup_nctoolbox

for it=1:length(tList)
    
    t=tList(it);  
   
    %Grib file name
    fname=sprintf('%s/%s/%s/namanl_218_%s_%s_000.grb',...
        'http://nomads.ncdc.noaa.gov/data/namanl',...
        datestr(t,'yyyymm'),...
        datestr(t,'yyyymmdd'),...
        datestr(t,'yyyymmdd'),...
        datestr(t,'HHMM')); 
    display(fname); 
    %fname='/home/server/student/homes/ipasmans/PhD/projects/boundary/NAM/namanl_218_20110401_0000_000.grb'; 
    
    %get overview content data set
    nco=ncgeodataset(fname); 
    
    %read Grib file
    val=nco{'Temperature_surface'}(1,grib.mask1,grib.mask2)-273.15; 
    val=squeeze(val); 
    grib.x=grib.x(grib.mask1,grib.mask2); 
    grib.y=grib.y(grib.mask1,grib.mask2); 
    
    
    %interpolate to lon output grid
    for i2=1:size(grib.x,2)
       valp1(:,i2)=interp1(grib.x(:,i2),val(:,i2),grd.x(:,1));  
       yp1(:,i2)=interp1(grib.x(:,i2),grib.y(:,i2),grd.x(:,1)); 
    end
    %interpolate to lat output grid
    for i1=1:size(valp1,1)
        valp2(i1,:)=interp1(yp1(i1,:),valp1(i1,:),grd.y(1,:)); 
    end            
    
    %plot
    close all; 
    figure();
    pcolor(grib.x,grib.y,double(val)); shading flat; 
    figure(); 
    pcolor(grd.x,grd.y,double(valp2)); shading flat; 
    
    %write to output
    ncwrite(out_file,varname,valp2,[1,1,it]);     
    
    
end