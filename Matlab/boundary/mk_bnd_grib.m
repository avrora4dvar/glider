%% Create surface forcing files from Grib files of the NAM-12km model
clear all; clc; 
cd('~/Matlab/boundary/'); 
addpath('~/Matlab/nctoolbox'); 
setup_nctoolbox

sig_file='~/Data/con/per_sig_bal.nc'; 
grid_file='~/Data/con/grd_ow2km_con.nc';
out_dir='~/Data/con/con2';
tList=[datenum('2011-02-28'):.25:datenum('2011-03-02')]; 
dateRef=datenum('2005-01-01'); 

%% Read grid files

%x:longitudes 0-360
%y:latitudes

%read grid coordinates output grid
grd.x=ncread(grid_file,'lon_rho'); grd.x=mod(grd.x,360); 
grd.y=ncread(grid_file,'lat_rho'); 
grd.mask=ncread(grid_file,'mask_rho'); 
grd.sig=ncread(sig_file,'temp',[1,1,40,1],[Inf,Inf,1,1]); grd.sig=grd.sig; 
gsize=[size(grd.x),length(tList)]; 

%read grid coordinate input grid
grib.grid=dlmread('NAM/latlon-g218.txt'); 
bsize=[max(grib.grid(:,1)),max(grib.grid(:,2))]; 
grib.x=mod(-reshape(grib.grid(:,4),bsize),360); 
grib.y=reshape(grib.grid(:,3),bsize); 

%wind-flow angle
grib.phi=nan(size(grib.x)); 
grib.phi(1:end-1,:)=atan( tan(deg2rad(grib.y(2:end,:)-grib.y(1:end-1,:)))./...
    sin(deg2rad(grib.x(2:end,:)-grib.x(1:end-1,:))) ); 

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

%make selection
grib.x=grib.x(grib.mask1,grib.mask2);
grib.y=grib.y(grib.mask1,grib.mask2);
grib.phi=grib.phi(grib.mask1,grib.mask2);


%% Create output files

timeunits='days since 2005-01-01'; 

%cloud coverage in [fraction 0-1]
fname=fullfile(out_dir,'frc_ow2km_Cloud.nc'); 
if exist(fname,'file'); delete(fname); end
nccreate(fname,'cloud_time','Format','classic','Dimensions',{'time',Inf}); 
nccreate(fname,'cloud','Format','classic','Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf},'Datatype','single'); 
ncwriteatt(fname,'cloud_time','units',timeunits);
ncwriteatt(fname,'cloud_time','long_name','cloud fraction time');
ncwriteatt(fname,'cloud','long_name','cloud fraction (0 to 1)');
ncwriteatt(fname,'cloud','units','nondimensional');
ncwriteatt(fname,'cloud','time','cloud_time');

%downward longwave radiation [W m-2]  (~200-300 W m-2)
fname=fullfile(out_dir,'frc_ow2km_Dlwrad.nc'); 
if exist(fname,'file'); delete(fname); end
nccreate(fname,'lrf_time','Format','classic','Dimensions',{'time',Inf}); 
nccreate(fname,'lwrad_down','Format','classic','Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf},'Datatype','single'); 
ncwriteatt(fname,'lrf_time','units',timeunits);
ncwriteatt(fname,'lrf_time','long_name','longwave radiation time');
ncwriteatt(fname,'lwrad_down','long_name','downward longwave radiation');
ncwriteatt(fname,'lwrad_down','units','W m-2');
ncwriteatt(fname,'lwrad_down','positive value','downward flux, heating');
ncwriteatt(fname,'lwrad_down','time','lrf_time');


%net longwave radiation [W m-2]
fname=fullfile(out_dir,'frc_ow2km_Lwrad.nc'); 
if exist(fname,'file'); delete(fname); end
nccreate(fname,'lrf_time','Format','classic','Dimensions',{'time',Inf}); 
nccreate(fname,'lwrad','Format','classic','Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf},'Datatype','single');
ncwriteatt(fname,'lrf_time','units',timeunits);
ncwriteatt(fname,'lrf_time','long_name','longwave radiation time');
ncwriteatt(fname,'lwrad','long_name','net longwave radiation');
ncwriteatt(fname,'lwrad','units','W m-2');
ncwriteatt(fname,'lwrad','positive value','downward flux, heating');
ncwriteatt(fname,'lwrad','time','lrf_time');


%net downward shortwave radiation [W m-2] (~100 W m-2)
fname=fullfile(out_dir,'frc_ow2km_Swrad.nc'); 
if exist(fname,'file'); delete(fname); end
nccreate(fname,'srf_time','Format','classic','Dimensions',{'time',Inf}); 
nccreate(fname,'swrad','Format','classic','Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf},'Datatype','single');
ncwriteatt(fname,'srf_time','units',timeunits);
ncwriteatt(fname,'srf_time','long_name','solar shortwave radiation time');
ncwriteatt(fname,'swrad','long_name','solar shortwave radiation');
ncwriteatt(fname,'swrad','units','W m-2');
ncwriteatt(fname,'swrad','time','srf_time');
ncwriteatt(fname,'swrad','positive value','downward flux, heating');
ncwriteatt(fname,'swrad','negative value','upward flux, cooling');


%air pressure [mb]
fname=fullfile(out_dir,'frc_ow2km_Pair.nc'); 
if exist(fname,'file'); delete(fname); end
nccreate(fname,'pair_time','Format','classic','Dimensions',{'time',Inf}); 
nccreate(fname,'Pair','Format','classic','Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf},'Datatype','single');
ncwriteatt(fname,'pair_time','units',timeunits);
ncwriteatt(fname,'pair_time','long_name','surface air pressure time');
ncwriteatt(fname,'Pair','long_name','surface air pressure');
ncwriteatt(fname,'Pair','units','milibar');
ncwriteatt(fname,'Pair','time','pair_time');

%relative humidity [% 0-100]
fname=fullfile(out_dir,'frc_ow2km_Qair.nc'); 
if exist(fname,'file'); delete(fname); end
nccreate(fname,'qair_time','Format','classic','Dimensions',{'time',Inf}); 
nccreate(fname,'Qair','Format','classic','Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf},'Datatype','single');
ncwriteatt(fname,'qair_time','units',timeunits);
ncwriteatt(fname,'qair_time','long_name','surface relative humidity time');
ncwriteatt(fname,'Qair','long_name','surface relative humidity');
ncwriteatt(fname,'Qair','units','milibar');
ncwriteatt(fname,'Qair','time','qair_time');

%air temperature [degC]
fname=fullfile(out_dir,'frc_ow2km_Tair.nc'); 
if exist(fname,'file'); delete(fname); end
nccreate(fname,'tair_time','Format','classic','Dimensions',{'time',Inf}); 
nccreate(fname,'Tair','Format','classic','Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf},'Datatype','single');
ncwriteatt(fname,'tair_time','units',timeunits);
ncwriteatt(fname,'tair_time','long_name','surface air temperature time');
ncwriteatt(fname,'Tair','long_name','surface air temperature');
ncwriteatt(fname,'Tair','units','Celsius');
ncwriteatt(fname,'Tair','time','tair_time');

%u-wind [m s-1]
fname=fullfile(out_dir,'frc_ow2km_wind.nc'); 
if exist(fname,'file'); delete(fname); end
nccreate(fname,'wind_time','Format','classic','Dimensions',{'time',Inf}); 
nccreate(fname,'Uwind','Format','classic','Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf},'Datatype','single');
nccreate(fname,'Vwind','Format','classic','Dimensions',{'xi_rho',gsize(1),'eta_rho',gsize(2),'time',Inf},'Datatype','single');
ncwriteatt(fname,'wind_time','units',timeunits);
ncwriteatt(fname,'wind_time','long_name','surface wind time');
ncwriteatt(fname,'Uwind','long_name','surface u-wind component');
ncwriteatt(fname,'Uwind','units','m s-1');
ncwriteatt(fname,'Uwind','time','wind_time');
ncwriteatt(fname,'Vwind','long_name','surface v-wind component');
ncwriteatt(fname,'Vwind','units','m s-1');
ncwriteatt(fname,'Vwind','time','wind_time');



%% Read fields

failList=[]; 
for it=1:length(tList)
   
    t=tList(it); 
    %Grib index file name
    iname=sprintf('%s/%s/%s/namanl_218_%s_%s_000.inv',...
        'https://nomads.ncdc.noaa.gov/data/namanl',...
        datestr(t,'yyyymm'),...
        datestr(t,'yyyymmdd'),...
        datestr(t,'yyyymmdd'),...
        datestr(t,'HHMM'));
    fname=sprintf('%s/%s/%s/namanl_218_%s_%s_000.grb',...
        'https://nomads.ncdc.noaa.gov/data/namanl',...
        datestr(t,'yyyymm'),...
        datestr(t,'yyyymmdd'),...
        datestr(t,'yyyymmdd'),...
        datestr(t,'HHMM'));
    display(fname);
    
    if exist('tmp.grb'); delete('tmp.grb'); end
    cd0=pwd(); 
    cd('~/Matlab/curl-master'); 
    command=sprintf('get_inv.pl %s | egrep "(%s)"|get_grib.pl %s %s',...
        iname,':TMP:2 m above gnd:|:RH:2 m above gnd:|:DSWRF:sfc:|:USWRF:sfc:|:DLWRF:sfc:|:ULWRF:sfc:|:PRMSL:MSL:|:TCDC:atmos col:|:UGRD:10 m above gnd:|:VGRD:10 m above gnd:',...
        fname,fullfile(cd0,'tmp.grib')); 
    unix(command); 
    cd(cd0); 
    
    if ~exist('tmp.grib','file')
        failList=[failList,it];
    else
        
        geo=ncgeodataset(fullfile(cd0,'tmp.grib'));
        
        %find grib points within model mask
        gribIn=interp2(grd.x',grd.y',double(grd.mask)',grib.x,grib.y); 
        gribIn(isnan(gribIn))=0; gribIn=gribIn(:)==1; 
        iAll=scatteredInterpolant(grib.x(:),grib.y(:),double(gribIn(:)),'linear'); 
        iIn=scatteredInterpolant(grib.x(gribIn),grib.y(gribIn),double(gribIn(gribIn)),'nearest'); 
        inGrd=iAll(grd.x,grd.y);
        inGrd=inGrd<=0.001|inGrd>=.999; 
        iAll.Method='natural'; 
                
        %temperature
        val0=geo{'Temperature_height_above_ground'}(1,1,grib.mask2,grib.mask1);
        val0=squeeze(double(val0))'-273.15; val1=nan(size(grd.mask)); 
        iAll.Values=val0(:); iIn.Values=val0(gribIn); 
        val1(inGrd)=iAll(grd.x(inGrd),grd.y(inGrd)); 
        val1(~inGrd)=iIn(grd.x(~inGrd),grd.y(~inGrd)); 
        ncwrite(fullfile(out_dir,'frc_ow2km_Tair.nc'),'Tair',val1,[1,1,it]);
        
        %relative-humidity
        val0=geo{'Relative_humidity_height_above_ground'}(1,1,grib.mask2,grib.mask1);
        val0=squeeze(double(val0))'; val1=nan(size(grd.mask)); 
        iAll.Values=val0(:); iIn.Values=val0(gribIn); 
        val1(inGrd)=iAll(grd.x(inGrd),grd.y(inGrd)); 
        val1(~inGrd)=iIn(grd.x(~inGrd),grd.y(~inGrd)); 
        ncwrite(fullfile(out_dir,'frc_ow2km_Qair.nc'),'Qair',val1,[1,1,it]); 
        
        %pressure
        val0=geo{'Pressure_reduced_to_MSL_msl'}(1,grib.mask2,grib.mask1);
        val0=squeeze(double(val0))'; val1=nan(size(grd.mask)); 
        iAll.Values=val0(:); iIn.Values=val0(gribIn); 
        val1(inGrd)=iAll(grd.x(inGrd),grd.y(inGrd)); 
        val1(~inGrd)=iIn(grd.x(~inGrd),grd.y(~inGrd)); 
        ncwrite(fullfile(out_dir,'frc_ow2km_Pair.nc'),'Pair',val1,[1,1,it]);
        
        %cloud cover
        val0=geo{'Total_cloud_cover_entire_atmosphere'}(1,grib.mask2,grib.mask1);
        val0=squeeze(double(val0))'/100; val1=nan(size(grd.mask)); 
        iAll.Values=val0(:); iIn.Values=val0(gribIn); 
        val1(inGrd)=iAll(grd.x(inGrd),grd.y(inGrd)); 
        val1(~inGrd)=iIn(grd.x(~inGrd),grd.y(~inGrd));
        ncwrite(fullfile(out_dir,'frc_ow2km_Cloud.nc'),'cloud',val1,[1,1,it]);
        
        %shortwave radiation
        val0d=geo{'downward_short_wave_rad_flux_surface'}(1,grib.mask2,grib.mask1);
        val0d=squeeze(double(val0d))';
        val0u=geo{'upward_short_wave_rad_flux_surface'}(1,grib.mask2,grib.mask1);
        val0u=squeeze(double(val0u))';
        val0=val0d-val0u; val1=nan(size(grd.mask)); 
        iAll.Values=val0(:); iIn.Values=val0(gribIn); 
        val1(inGrd)=iAll(grd.x(inGrd),grd.y(inGrd)); 
        val1(~inGrd)=iIn(grd.x(~inGrd),grd.y(~inGrd));
        ncwrite(fullfile(out_dir,'frc_ow2km_Swrad.nc'),'swrad',val1,[1,1,it]);
        
        %down longwave
        val0d=geo{'downward_long_wave_rad_flux_surface'}(1,grib.mask2,grib.mask1);
        val0d=squeeze(double(val0d))';
        val0=val0d; val1=nan(size(grd.mask)); 
        iAll.Values=val0(:); iIn.Values=val0(gribIn); 
        val1(inGrd)=iAll(grd.x(inGrd),grd.y(inGrd)); 
        val1(~inGrd)=iIn(grd.x(~inGrd),grd.y(~inGrd));
        ncwrite(fullfile(out_dir,'frc_ow2km_Dlwrad.nc'),'lwrad_down',val1,[1,1,it]);
        
        %net longwave
        val0d=geo{'downward_long_wave_rad_flux_surface'}(1,grib.mask2,grib.mask1);
        val0d=squeeze(double(val0d))';
        val0u=geo{'upward_long_wave_rad_flux_surface'}(1,grib.mask2,grib.mask1);
        val0u=squeeze(double(val0u))';
        val0=val0d-val0u; val1=nan(size(grd.mask)); 
        iAll.Values=val0(:); iIn.Values=val0(gribIn); 
        val1(inGrd)=iAll(grd.x(inGrd),grd.y(inGrd)); 
        val1(~inGrd)=iIn(grd.x(~inGrd),grd.y(~inGrd));
        ncwrite(fullfile(out_dir,'frc_ow2km_Lwrad.nc'),'lwrad',val1,[1,1,it])
        
        %calculate-wind directions
        val0uL=geo{'u-component_of_wind_height_above_ground'}(1,1,grib.mask2,grib.mask1);
        val0uL=squeeze(double(val0uL))';
        val0vL=geo{'v-component_of_wind_height_above_ground'}(1,1,grib.mask2,grib.mask1);
        val0vL=squeeze(double(val0vL))';
        
        val0u=cos(grib.phi).*val0uL-sin(grib.phi).*val0vL;
        val0v=cos(grib.phi).*val0vL+sin(grib.phi).*val0uL;
        
        val1=nan(size(grd.mask)); 
        iAll.Values=val0u(:); iIn.Values=val0u(gribIn); 
        val1(inGrd)=iAll(grd.x(inGrd),grd.y(inGrd)); 
        val1(~inGrd)=iIn(grd.x(~inGrd),grd.y(~inGrd));
        ncwrite(fullfile(out_dir,'frc_ow2km_wind.nc'),'Uwind',val1,[1,1,it]);
        
        val1=nan(size(grd.mask)); 
        iAll.Values=val0v(:); iIn.Values=val0v(gribIn); 
        val1(inGrd)=iAll(grd.x(inGrd),grd.y(inGrd)); 
        val1(~inGrd)=iIn(grd.x(~inGrd),grd.y(~inGrd));
        ncwrite(fullfile(out_dir,'frc_ow2km_wind.nc'),'Vwind',val1,[1,1,it]);
        
    end %if exist

end

error('ok'); 

%% interpolate for missing values

itList=[1:length(tList)]; 
for it=1:length(tList)
   if ismember(it,failList)
       %find at same moment of days
       ib1=find(itList<it & ~ismember(itList,failList) & mod(tList-tList(it),1)==0,1,'last'); 
       ib2=find(itList>it & ~ismember(itList,failList) & mod(tList-tList(it),1)==0,1,'first');
       w=(tList(it)-tList(ib1))/(tList(ib2)-tList(ib1)); 
       
       val0=ncread(fullfile(out_dir,'frc_ow2km_Tair.nc'),'Tair',[1,1,ib1],[Inf,Inf,1]);
       val1=ncread(fullfile(out_dir,'frc_ow2km_Tair.nc'),'Tair',[1,1,ib2],[Inf,Inf,1]);
       ncwrite(fullfile(out_dir,'frc_ow2km_Tair.nc'),'Tair',(1-w)*val0+w*val1,[1,1,it]);
       
       val0=ncread(fullfile(out_dir,'frc_ow2km_Qair.nc'),'Qair',[1,1,ib1],[Inf,Inf,1]);
       val1=ncread(fullfile(out_dir,'frc_ow2km_Qair.nc'),'Qair',[1,1,ib2],[Inf,Inf,1]);
       ncwrite(fullfile(out_dir,'frc_ow2km_Qair.nc'),'Qair',(1-w)*val0+w*val1,[1,1,it]);       

       
       val0=ncread(fullfile(out_dir,'frc_ow2km_Dlwrad.nc'),'lwrad_down',[1,1,ib1],[Inf,Inf,1]);
       val1=ncread(fullfile(out_dir,'frc_ow2km_Dlwrad.nc'),'lwrad_down',[1,1,ib2],[Inf,Inf,1]);
       ncwrite(fullfile(out_dir,'frc_ow2km_Dlwrad.nc'),'lwrad_down',(1-w)*val0+w*val1,[1,1,it]);
       
       val0=ncread(fullfile(out_dir,'frc_ow2km_Lwrad.nc'),'lwrad',[1,1,ib1],[Inf,Inf,1]);
       val1=ncread(fullfile(out_dir,'frc_ow2km_Lwrad.nc'),'lwrad',[1,1,ib2],[Inf,Inf,1]);
       ncwrite(fullfile(out_dir,'frc_ow2km_Lwrad.nc'),'lwrad',(1-w)*val0+w*val1,[1,1,it]);
       
       val0=ncread(fullfile(out_dir,'frc_ow2km_Swrad.nc'),'swrad',[1,1,ib1],[Inf,Inf,1]);
       val1=ncread(fullfile(out_dir,'frc_ow2km_Swrad.nc'),'swrad',[1,1,ib2],[Inf,Inf,1]);
       ncwrite(fullfile(out_dir,'frc_ow2km_Swrad.nc'),'swrad',(1-w)*val0+w*val1,[1,1,it]);
       
       %find closest in time
       ib1=find(itList<it & ~ismember(itList,failList) ,1,'last'); 
       ib2=find(itList>it & ~ismember(itList,failList) ,1,'first');
       w=(tList(it)-tList(ib1))/(tList(ib2)-tList(ib1)); 
       
       val0=ncread(fullfile(out_dir,'frc_ow2km_Pair.nc'),'Pair',[1,1,ib1],[Inf,Inf,1]);
       val1=ncread(fullfile(out_dir,'frc_ow2km_Pair.nc'),'Pair',[1,1,ib2],[Inf,Inf,1]);
       ncwrite(fullfile(out_dir,'frc_ow2km_Pair.nc'),'Pair',(1-w)*val0+w*val1,[1,1,it]);
       
       val0=ncread(fullfile(out_dir,'frc_ow2km_Cloud.nc'),'cloud',[1,1,ib1],[Inf,Inf,1]);
       val1=ncread(fullfile(out_dir,'frc_ow2km_Cloud.nc'),'cloud',[1,1,ib2],[Inf,Inf,1]);
       ncwrite(fullfile(out_dir,'frc_ow2km_Cloud.nc'),'cloud',(1-w)*val0+w*val1,[1,1,it]);
       
       val0=ncread(fullfile(out_dir,'frc_ow2km_wind.nc'),'Uwind',[1,1,ib1],[Inf,Inf,1]);
       val1=ncread(fullfile(out_dir,'frc_ow2km_wind.nc'),'Uwind',[1,1,ib2],[Inf,Inf,1]);
       ncwrite(fullfile(out_dir,'frc_ow2km_wind.nc'),'Uwind',(1-w)*val0+w*val1,[1,1,it]);
       
       val0=ncread(fullfile(out_dir,'frc_ow2km_wind.nc'),'Vwind',[1,1,ib1],[Inf,Inf,1]);
       val1=ncread(fullfile(out_dir,'frc_ow2km_wind.nc'),'Vwind',[1,1,ib2],[Inf,Inf,1]);
       ncwrite(fullfile(out_dir,'frc_ow2km_wind.nc'),'Vwind',(1-w)*val0+w*val1,[1,1,it]);
       
   end
end

%% test for NaN; 

val=ncread(fullfile(out_dir,'frc_ow2km_wind.nc'),'Uwind'); 
nanval.uwind=find( any(reshape(isnan(val),[],size(val,3)),1) ); 
val=ncread(fullfile(out_dir,'frc_ow2km_wind.nc'),'Vwind'); 
nanval.vwind=find( any(reshape(isnan(val),[],size(val,3)),1) ); 
val=ncread(fullfile(out_dir,'frc_ow2km_Tair.nc'),'Tair'); 
nanval.tair=find( any(reshape(isnan(val),[],size(val,3)),1) ); 
val=ncread(fullfile(out_dir,'frc_ow2km_Qair.nc'),'Qair'); 
nanval.qair=find( any(reshape(isnan(val),[],size(val,3)),1) ); 
val=ncread(fullfile(out_dir,'frc_ow2km_Pair.nc'),'Pair'); 
nanval.pair=find( any(reshape(isnan(val),[],size(val,3)),1) ); 
val=ncread(fullfile(out_dir,'frc_ow2km_Dlwrad.nc'),'lwrad_down'); 
nanval.dlwrad=find( any(reshape(isnan(val),[],size(val,3)),1) ); 
val=ncread(fullfile(out_dir,'frc_ow2km_Lwrad.nc'),'lwrad');
nanval.lwrad=find( any(reshape(isnan(val),[],size(val,3)),1) ); 
val=ncread(fullfile(out_dir,'frc_ow2km_Swrad.nc'),'swrad'); 
nanval.swrad=find( any(reshape(isnan(val),[],size(val,3)),1) ); 
val=ncread(fullfile(out_dir,'frc_ow2km_Cloud.nc'),'cloud'); 
nanval.cloud=find( any(reshape(isnan(val),[],size(val,3)),1) ); 


%% times

ncwrite(fullfile(out_dir,'frc_ow2km_wind.nc'),'wind_time',tList-dateRef); 
ncwrite(fullfile(out_dir,'frc_ow2km_Cloud.nc'),'cloud_time',tList-dateRef); 
ncwrite(fullfile(out_dir,'frc_ow2km_Dlwrad.nc'),'lrf_time',tList-dateRef); 
ncwrite(fullfile(out_dir,'frc_ow2km_Lwrad.nc'),'lrf_time',tList-dateRef);
ncwrite(fullfile(out_dir,'frc_ow2km_Swrad.nc'),'srf_time',tList-dateRef); 
ncwrite(fullfile(out_dir,'frc_ow2km_Tair.nc'),'tair_time',tList-dateRef); 
ncwrite(fullfile(out_dir,'frc_ow2km_Qair.nc'),'qair_time',tList-dateRef); 
ncwrite(fullfile(out_dir,'frc_ow2km_Pair.nc'),'pair_time',tList-dateRef); 



