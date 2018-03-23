% Lana ctz implementation
% 6/23/04: fills in file frc_bulk.nc 
% (filling into a template for frc_bulk.nc prepared with the CDL file)
% load

% Fields:
% Uwind, Vwind, Tair, Qair (RH), swrad, cloud + associated times
% ctz_ts
%
function [status]=heat_primes_ctz(reini);
%
wrkdir='/home/server/student/homes/ipasmans/PhD/projects/boundary/';
ldir=[wrkdir 'matlib/'];
path(ldir,path);
path([ldir 'AIR_SEA'],path);
path([ldir 'FUNCTIONS'],path);
%
date_ref='01-Jan-2005'; % ref date and time in roms application
day_ref=0;
timeunits=['days since ' datestr(date_ref,'yyyy-mm-dd HH:MM:SS')];
%
dayn=floor(now); %REMOVE -1!
today=datestr(dayn);
year=str2num(today(8:11));
day_now=floor(now-datenum(year,1,1));
%
Gname=[wrkdir 'Prm/grd_ow2km_03_smooth04_new.nc'];
mask_rho=ncread(Gname,'mask_rho');
lon=ncread(Gname,'lon_rho');
lat=ncread(Gname,'lat_rho');
%
[Mr,Lr]=size(mask_rho); % same as [eta_rho,xi_rho]
eta_rho=Lr;xi_rho=Mr;
%
if reini==1,
 load CTZ_ts_var0;
 fnameC=[ wrkdir 'Prm/frc_ow2km_Cloud0.nc'];
 fnameP=[ wrkdir 'Prm/frc_ow2km_Pair0.nc'];
 fnameQ=[ wrkdir 'Prm/frc_ow2km_Qair0.nc'];
 fnameS=[ wrkdir 'Prm/frc_ow2km_Swrad0.nc'];
 fnameT=[ wrkdir 'Prm/frc_ow2km_Tair0.nc'];
 fprintf('REINITIALIZATION CASE, stage 1\n');
elseif reini==2,
 load ctz_ts_1;
 fnameC=[ wrkdir 'Prm/frc_ow2km_Cloud_1.nc'];
 fnameP=[ wrkdir 'Prm/frc_ow2km_Pair_1.nc'];
 fnameQ=[ wrkdir 'Prm/frc_ow2km_Qair_1.nc'];
 fnameS=[ wrkdir 'Prm/frc_ow2km_Swrad_1.nc'];
 fnameT=[ wrkdir 'Prm/frc_ow2km_Tair_1.nc'];
 fprintf('REINITIALIZATION CASE, stage 2\n');
else
 load CTZ_ts_var;
 fnameC=[ wrkdir 'Prm/frc_ow2km_Cloud.nc'];
 fnameP=[ wrkdir 'Prm/frc_ow2km_Pair.nc'];
 fnameQ=[ wrkdir 'Prm/frc_ow2km_Qair.nc'];
 fnameS=[ wrkdir 'Prm/frc_ow2km_Swrad.nc'];
 fnameT=[ wrkdir 'Prm/frc_ow2km_Tair.nc'];
 fprintf('REGULAR RUN\n');
end
%
% Initialize the forcing files: (2014)
status=define_roms_cloud_file(fnameC,xi_rho,eta_rho,Gname,timeunits);
status=define_roms_pair_file(fnameP,xi_rho,eta_rho,Gname,timeunits);
status=define_roms_qair_file(fnameQ,xi_rho,eta_rho,Gname,timeunits);
status=define_roms_swrad_file(fnameS,xi_rho,eta_rho,Gname,timeunits);
status=define_roms_tair_file(fnameT,xi_rho,eta_rho,Gname,timeunits);
%
T1=time(1);
T2=time(2); %T2=time(end);
t0=time(1:2); %t0=time;
%t0=timef;
nt=length(t0);
ea=1; % <- output at times t0(1:ea:end);
t0loc=t0;     
yd=[T1:T2]';
%
zr=10;   % wind mes. sensor height
zt=2.;   % air temp. sensor height
zq=2.;   % Rel. Hum. sensor height
%
mlon1=min(min(lon))-0.1;
mlon2=max(max(lon))+0.1;
mlat1=min(min(lat))-0.1;
mlat2=max(max(lat))+0.1;
[dum,k1]=min(abs(mlon-mlon1)+abs(mlat-mlat1));
[dum,k2]=min(abs(mlon-mlon2)+abs(mlat-mlat1));
[dum,k3]=min(abs(mlon-mlon1)+abs(mlat-mlat2));
[dum,k4]=min(abs(mlon-mlon2)+abs(mlat-mlat2));
%
% append with corners 
mlon=[mlon;mlon1;mlon2;mlon1;mlon2];
mlat=[mlat;mlat1;mlat1;mlat2;mlat2];
%%%%%%%%%%%%%%%%%%%%%
% air temp ETA
%%%%%%%%%%%%%%%%%%%%%
%[Ta0,tTalp]=OSUlpAK(T_2m_ts,time);
Ta0=T_2m_ts;tTalp=time;

% append with closest to corners values
Ta0=[Ta0;Ta0(k1,:);Ta0(k2,:);Ta0(k3,:);Ta0(k4,:)];
% interpolate on lon,lat
k=1;
for it=1:ea:nt
 fprintf('Interpolating Air Temperature on orw2km grid, step %d of %d...',it,nt);
 Tair_nc(:,:,k)=griddata(mlon,mlat,Ta0(:,it),lon,lat);k=k+1;
 fprintf('done\n');
end
%%%%%%%%%%%%%%%%%%%%%%
%% sea surface pressure
Pa=0.01*SurfPr_ts;              % <- no low passed
%[Pa0,tPalp]=OSUlpAK(Pa,time);
Pa0=Pa;tPalp=time;
Pa0=[Pa0;Pa0(k1,:);Pa0(k2,:);Pa0(k3,:);Pa0(k4,:)];
% interpolate on lon,lat
k=1;
for it=1:ea:nt
 fprintf('Interpolating Air Pressure on orw2km grid, step %d of %d...',it,nt);
 Pair_nc(:,:,k)=griddata(mlon,mlat,Pa0(:,it),lon,lat);k=k+1;
 fprintf('done\n');
end
%%%%%%%%%%%%%%%%%%%%%
% Rel. humidity data:
%%%%%%%%%%%%%%%%%%%%%
RH=RelHum_ts;RH(find(RH>100))=100;
%[RH0,tRHlp]=OSUlpAK(RH,time);
RH0=RH;tRHlp=time;
RH0=[RH0;RH0(k1,:);RH0(k2,:);RH0(k3,:);RH0(k4,:)];
% interpolate on lon,lat
k=1;
for it=1:ea:nt
 fprintf('Interpolating Rel Humidity on orw2km grid, step %d of %d...',it,nt);
 Qair_nc(:,:,k)=griddata(mlon,mlat,RH0(:,it),lon,lat);k=k+1;
 fprintf('done\n');
end
%%%%%%%%%%%%%%%%%%%%%
% Short wave radiation:
%%%%%%%%%%%%%%%%%%%%%
Qsw=DownSwave_ts;
%[Qsw0,tswlp]=OSUlpAK(Qsw,time);
Qsw0=Qsw;tswlp=time;
Qsw0=[Qsw0;Qsw0(k1,:);Qsw0(k2,:);Qsw0(k3,:);Qsw0(k4,:)];
% interpolate on lon,lat
k=1;
for it=1:ea:nt
 fprintf('Interpolating short wave radiation on orw2km grid, step %d of %d...',it,nt);
 swrad_nc(:,:,k)=griddata(mlon,mlat,Qsw0(:,it),lon,lat);k=k+1;
 fprintf('done\n');
end
%%%%%%%%%%%%%%%%%%%%%
% Cloud cover:
%%%%%%%%%%%%%%%%%%%%%
C0loc=0.01*CloudCov_ts;
C0loc(find(C0loc<0))=0;
C0loc(find(C0loc>1))=1;
%[C0,tlp]=OSUlpAK(C0loc,time);
C0=C0loc;tlp=time;
C0=[C0;C0(k1,:);C0(k2,:);C0(k3,:);C0(k4,:)];
% interpolate on lon,lat
k=1;
for it=1:ea:nt
 fprintf('Interpolating Cloud Cover on orw2km grid, step %d of %d...',it,nt);
 cloud_nc(:,:,k)=griddata(mlon,mlat,C0(:,it),lon,lat);k=k+1;
 fprintf('done\n');
end
%
t_nc=t0(1:ea:end); 
nt=length(t_nc);
ncwrite(fnameT,'tair_time', t_nc);
ncwrite(fnameQ,'qair_time', t_nc);
ncwrite(fnameP,'pair_time', t_nc);
ncwrite(fnameS,'srf_time',  t_nc);
ncwrite(fnameC,'cloud_time',t_nc);
for rec=1:nt
 tmp=Tair_nc(:,:,rec);
 ncwrite(fnameT,'Tair',tmp,[1,1,rec]);
 tmp=Qair_nc(:,:,rec);
 ncwrite(fnameQ,'Qair',tmp,[1,1,rec]);
 tmp=Pair_nc(:,:,rec);
 ncwrite(fnameP,'Pair',tmp,[1,1,rec]);
 tmp=swrad_nc(:,:,rec);
 ncwrite(fnameS,'swrad',tmp,[1,1,rec]);
 tmp=cloud_nc(:,:,rec);
 ncwrite(fnameC,'cloud',tmp,[1,1,rec]);
end

