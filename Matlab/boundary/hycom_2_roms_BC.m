% version 2: read HYCOM fields from multiple files:

clear all;
romsgrid='/home/server/student/homes/ipasmans/Matlab/projects/smooth_bathy/grd_ow2km_con.nc';
if exist('CTZ_ts_var.mat','file'),
  load CTZ_ts_var;
else
  load CTZ_ts_var_y;
end
date_ref='01-Jan-2005';
timec=[datenum('2011-01-01'):datenum('2012-01-01')]; 
time=timec-datenum(date_ref); 
tnow=timec(1); tend=timec(end); 
roms_dates={datestr(floor(timec(1))-1),datestr(tend)};
%%%%%roms_dates={'12-Mar-2015','31-Mar-2015'};
%
 N=40;
 Vtransform=2;
 Vstretching=4;
 theta_s=8;     % parameter for stretching near surface
 theta_b=3;   % parameter for stretching near bottom
 Tcline=50;     % thermocline depth
%
 day_ref=0;
%
 fnameout='/home/server/student/homes/ipasmans/Data/bnd_ow2km_con_hycom_2011.nc'; %%% !!!!
%
 new_file=1;  % if 0, continue write to the old file (e.g., when processing was interrupted)
 it_start=1; % new_file=1, this will be reset to 1
  
 north_yes=1;
 south_yes=1;
 west_yes=1;
 east_yes=0;


if new_file==0
 reply=input('Are you sure you want to append to the existing file? Y/N [Y]:','s');
 if isempty(reply)
  reply = 'Y';
 end
 if reply=='N' or reply=='n'
  return;
 end
end

if new_file
 it_start=1;
end

% Hycom URL address:
% HYCOM + NCODA Global 1/12° Analysis (GLBa0.08) 	Print	
%    91.1 (Apr-2014 to Present)
%    91.0 (Aug-2013 to Apr-2014)
%    90.9 (Jan-2011 to Aug-2013)
%    90.8 (May-2009 to Jan-2011)
%    90.6 (Sep-2008 to May-2009)

urlhead='http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_';
urlhead1='http://tds.hycom.org/thredds/dodsC/GLBa0.08/latest/';

urlends={'90.6/2008/';...
         '90.6/2009/';...
         '90.8/2009/';...
         '90.8/2010/';...
         '90.8/2011/';...
         '90.9/2011/';...
         '90.9/2012/';...
         '90.9/2013/';...
         '91.0/2013/';...
         '91.0/2014/';...
         '91.1/2014/';...
         '91.1/2015/'};

nf=length(urlends);

% Find the url dir with the first record

roms_day_1=datenum(roms_dates{1});
roms_day_2=datenum(roms_dates{2});
roms_days=[roms_day_1:roms_day_2];
roms_day_start=roms_days(it_start);
nt=length(roms_days);

kf1=0;
kf=1; % Lana's fix, 2014 Dec 28
while kf1==0
 if kf<=nf
  urldir=[urlhead urlends{kf}]; % use the first file for mask definition etc.
  url_2d1=[urldir '2d']
  hycom_date=ncread(url_2d1,'Date');
  day_h=datenum(int2str(hycom_date),'yyyymmdd');
  it1=find(day_h==roms_day_start);
  if ~isempty(it1)
   disp(['The first record: ' int2str(it1) ' in ' urlends{kf}]);
   kf1=kf;
  end
  kf=kf+1;
 else
  error('the first roms date is not within the range of HYCOM');
 end
end

urldir=[urlhead urlends{kf1}];

% Find the url dir with the last record


kf2=0;
kf=1;
while kf2==0
 if kf<=nf
  urldir=[urlhead urlends{kf}]; % use the first file for mask definition etc.
  url_2d2=[urldir '2d'];
  hycom_date=ncread(url_2d2,'Date');
  day_h=datenum(int2str(hycom_date),'yyyymmdd');
  it2=find(day_h==roms_day_2);
  if ~isempty(it2)
   disp(['The last record: ' int2str(it2) ' in ' urlends{kf}]);
   kf2=kf;
  end
  kf=kf+1;
 else
%  urldir=urlhead1;
%  url_2d2=[urldir '2d'];
%  hycom_date=ncread(url_2d2,'Date');
%  day_h=datenum(int2str(hycom_date),'yyyymmdd');
%  it2=find(day_h==roms_day_2);
%  if ~isempty(it2)
%   disp(['The last record: ' int2str(it2) ' in latest']);
%   kf2=kf+1;
%  else
   error('the last roms date is not within the range of HYCOM');
%  end
 end
end

% roms coordinates, depth
lon_rho=ncread(romsgrid,'lon_rho');
lat_rho=ncread(romsgrid,'lat_rho');
lon_u=ncread(romsgrid,'lon_u');
lat_u=ncread(romsgrid,'lat_u');
lon_v=ncread(romsgrid,'lon_v');
lat_v=ncread(romsgrid,'lat_v');
h=ncread(romsgrid,'h');

[xi_rho,eta_rho]=size(lon_rho);

% Define the output file:
timeunits=['days since ' datestr(date_ref,'yyyy-mm-dd HH:MM:SS')];

if new_file
 tmp=define_roms_bc_file(fnameout,xi_rho,eta_rho,N,...
    theta_s,theta_b,Tcline,Vtransform,Vstretching,romsgrid,timeunits,...
    north_yes,south_yes,west_yes,east_yes);
end

z3d=zeros(xi_rho,eta_rho,N);
zeta0=zeros(xi_rho,eta_rho);
kgrid=0; % 0 for rho, 1 for W points
column=1;
plt=0;
for j=1:eta_rho
 [z0j,sc,Cs]=scoord_new(h,zeta0,theta_s,theta_b,Tcline,N,...
                        kgrid,column,j,Vtransform,Vstretching,plt);
 z3d(:,j,:)=reshape(z0j,[xi_rho 1 N]);
end

z_w=zeros(xi_rho,eta_rho,N+1);
zeta0=zeros(xi_rho,eta_rho);
kgrid=1; % 0 for rho, 1 for W points
column=1;
plt=0;
for j=1:eta_rho
 [z0j,sc,Cs]=scoord_new(h,zeta0,theta_s,theta_b,Tcline,N,...
                        kgrid,column,j,Vtransform,Vstretching,plt);
 z_w(:,j,:)=reshape(z0j,[xi_rho 1 N+1]);
end
Hz=z_w(:,:,2:end)-z_w(:,:,1:end-1); % depth of each layer (3D)

% HYCOM coordinates (read from file):
disp('load HYCOM coordinates (from a local file)...');
load hycom_lon_lat X Y Lon Lat;
Lon=double(Lon);
Lat=double(Lat);
X=double(X);
Y=double(Y);

url_uvel=[urldir 'uvel'];
z=-ncread(url_uvel,'Depth');
z=double(z);
z=flipud(z);
z=[-14000;z];
nz=length(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Northern boundary, grid etc.:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if north_yes
 disp('Northern Boundary...');
 lon_r_n=lon_rho(:,end);
 lat_r_n=lat_rho(:,end);
 lon_u_n=lon_u(:,end);
 lat_u_n=lat_u(:,end);
 lon_v_n=lon_v(:,end);
 lat_v_n=lat_v(:,end);

 hr_n=h(:,end);
 hu_n=0.5*(h(1:end-1,end)+h(2:end,end));
 hv_n=mean(h(:,end-1:end),2);

 zr_n=z3d(:,end,:);
 zu_n=0.5*(z3d(1:end-1,end,:)+z3d(2:end,end,:));
 zv_n=0.5*(z3d(:,end-1,:)+z3d(:,end,:));
 zr_n=squeeze(zr_n);
 zu_n=squeeze(zu_n);
 zv_n=squeeze(zv_n);

 Hz_u_n=0.5*(Hz(1:end-1,end,:)+Hz(2:end,end,:));
 Hz_v_n=0.5*(Hz(:,end-1,:)+Hz(:,end,:));
 Hz_u_n=squeeze(Hz_u_n);
 Hz_v_n=squeeze(Hz_v_n);

 disp('Find HYCOM subgrid bounds...');
 [i1_n,i2_n,j1_n,j2_n,xbr_n,ybr_n,XX_n,YY_n]=...
     find_hycom_bounds(lon_r_n,lat_r_n,X,Y,Lon,Lat);
 [i1_n,i2_n,j1_n,j2_n,xbu_n,ybu_n,XX_n,YY_n]=...
     find_hycom_bounds(lon_u_n,lat_u_n,X,Y,Lon,Lat);
 [i1_n,i2_n,j1_n,j2_n,xbv_n,ybv_n,XX_n,YY_n]=...
     find_hycom_bounds(lon_v_n,lat_v_n,X,Y,Lon,Lat);
 ii_n=[i1_n:i2_n];
 jj_n=[j1_n:j2_n];
 nx_n=length(ii_n);
 ny_n=length(jj_n);

 xbr2_n=repmat(xbr_n,[1 N]);
 ybr2_n=repmat(ybr_n,[1 N]);
 xbu2_n=repmat(xbu_n,[1 N]);
 ybu2_n=repmat(ybu_n,[1 N]);
 xbv2_n=repmat(xbv_n,[1 N]);
 ybv2_n=repmat(ybv_n,[1 N]);

end % north_yes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Southern boundary, grid etc.:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Southern Boundary...');
if south_yes
 lon_r_s=lon_rho(:,1);
 lat_r_s=lat_rho(:,1);
 lon_u_s=lon_u(:,1);
 lat_u_s=lat_u(:,1);
 lon_v_s=lon_v(:,1);
 lat_v_s=lat_v(:,1);

 hr_s=h(:,1);
 hu_s=0.5*(h(1:end-1,1)+h(2:end,1));
 hv_s=mean(h(:,1:2),2);

 zr_s=z3d(:,1,:);
 zu_s=0.5*(z3d(1:end-1,1,:)+z3d(2:end,1,:));
 zv_s=0.5*(z3d(:,1,:)+z3d(:,2,:));
 zr_s=squeeze(zr_s);
 zu_s=squeeze(zu_s);
 zv_s=squeeze(zv_s);

 Hz_u_s=0.5*(Hz(1:end-1,1,:)+Hz(2:end,1,:));
 Hz_v_s=0.5*(Hz(:,1,:)+Hz(:,2,:));
 Hz_u_s=squeeze(Hz_u_s);
 Hz_v_s=squeeze(Hz_v_s);

 disp('Find HYCOM subgrid bounds...');
 [i1_s,i2_s,j1_s,j2_s,xbr_s,ybr_s,XX_s,YY_s]=...
     find_hycom_bounds(lon_r_s,lat_r_s,X,Y,Lon,Lat);
 [i1_s,i2_s,j1_s,j2_s,xbu_s,ybu_s,XX_s,YY_s]=...
     find_hycom_bounds(lon_u_s,lat_u_s,X,Y,Lon,Lat);
 [i1_s,i2_s,j1_s,j2_s,xbv_s,ybv_s,XX_s,YY_s]=...
     find_hycom_bounds(lon_v_s,lat_v_s,X,Y,Lon,Lat);
 ii_s=[i1_s:i2_s];
 jj_s=[j1_s:j2_s];
 nx_s=length(ii_s);
 ny_s=length(jj_s);

 xbr2_s=repmat(xbr_s,[1 N]);
 ybr2_s=repmat(ybr_s,[1 N]);
 xbu2_s=repmat(xbu_s,[1 N]);
 ybu2_s=repmat(ybu_s,[1 N]);
 xbv2_s=repmat(xbv_s,[1 N]);
 ybv2_s=repmat(ybv_s,[1 N]);

end % south_yes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Western boundary, grid etc.:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Western Boundary...');

if west_yes
 lon_r_w=lon_rho(1,:);
 lat_r_w=lat_rho(1,:);
 lon_u_w=lon_u(1,:);
 lat_u_w=lat_u(1,:);
 lon_v_w=lon_v(1,:);
 lat_v_w=lat_v(1,:);

 hr_w=h(1,:);
 hu_w=0.5*(h(1,:)+h(2,:));
 hv_w=0.5*(h(1,1:end-1)+h(1,2:end));

 zr_w=z3d(1,:,:);
 zu_w=0.5*(z3d(1,:,:)+z3d(2,:,:));
 zv_w=0.5*(z3d(1,1:end-1,:)+z3d(1,2:end,:));
 zr_w=squeeze(zr_w);
 zu_w=squeeze(zu_w);
 zv_w=squeeze(zv_w);

 Hz_u_w=0.5*(Hz(1,:,:)+Hz(2,:,:));
 Hz_v_w=0.5*(Hz(1,1:end-1,:)+Hz(1,2:end,:));
 Hz_u_w=squeeze(Hz_u_w);
 Hz_v_w=squeeze(Hz_v_w);

 disp('Find HYCOM subgrid bounds...');
 [i1_w,i2_w,j1_w,j2_w,xbr_w,ybr_w,XX_w,YY_w]=...
     find_hycom_bounds(lon_r_w,lat_r_w,X,Y,Lon,Lat);
 [i1_w,i2_w,j1_w,j2_w,xbu_w,ybu_w,XX_w,YY_w]=...
     find_hycom_bounds(lon_u_w,lat_u_w,X,Y,Lon,Lat);
 [i1_w,i2_w,j1_w,j2_w,xbv_w,ybv_w,XX_w,YY_w]=...
     find_hycom_bounds(lon_v_w,lat_v_w,X,Y,Lon,Lat);
 ii_w=[i1_w:i2_w];
 jj_w=[j1_w:j2_w];
 nx_w=length(ii_w);
 ny_w=length(jj_w);

 xbr2_w=repmat(xbr_w',[1 N]);
 ybr2_w=repmat(ybr_w',[1 N]);
 xbu2_w=repmat(xbu_w',[1 N]);
 ybu2_w=repmat(ybu_w',[1 N]);
 xbv2_w=repmat(xbv_w',[1 N]);
 ybv2_w=repmat(ybv_w',[1 N]);

end % west_yes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eastern boundary, grid etc.:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if east_yes

 disp('Eastern Boundary...');
 lon_r_e=lon_rho(end,:);
 lat_r_e=lat_rho(end,:);
 lon_u_e=lon_u(end,:);
 lat_u_e=lat_u(end,:);
 lon_v_e=lon_v(end,:);
 lat_v_e=lat_v(end,:);

 hr_e=h(end,:);
 hu_e=0.5*(h(end-1,:)+h(end,:));
 hv_e=0.5*(h(end,1:end-1)+h(end,2:end));

 zr_e=z3d(end,:,:);
 zu_e=0.5*(z3d(end-1,:,:)+z3d(end,:,:));
 zv_e=0.5*(z3d(end,1:end-1,:)+z3d(end,2:end,:));
 zr_e=squeeze(zr_e);
 zu_e=squeeze(zu_e);
 zv_e=squeeze(zv_e);

 Hz_u_e=0.5*(Hz(end-1,:,:)+Hz(end,:,:));
 Hz_v_e=0.5*(Hz(end,1:end-1,:)+Hz(end,2:end,:));
 Hz_u_e=squeeze(Hz_u_e);
 Hz_v_e=squeeze(Hz_v_e);

 disp('Find HYCOM subgrid bounds...');
 [i1_e,i2_e,j1_e,j2_e,xbr_e,ybr_e,XX_e,YY_e]=...
     find_hycom_bounds(lon_r_e,lat_r_e,X,Y,Lon,Lat);
 [i1_e,i2_e,j1_e,j2_e,xbu_e,ybu_e,XX_e,YY_e]=...
     find_hycom_bounds(lon_u_e,lat_u_e,X,Y,Lon,Lat);
 [i1_e,i2_e,j1_e,j2_e,xbv_e,ybv_e,XX_e,YY_e]=...
     find_hycom_bounds(lon_v_e,lat_v_e,X,Y,Lon,Lat);
 ii_e=[i1_e:i2_e];
 jj_e=[j1_e:j2_e];
 nx_e=length(ii_e);
 ny_e=length(jj_e);

 xbr2_e=repmat(xbr_e',[1 N]);
 ybr2_e=repmat(ybr_e',[1 N]);
 xbu2_e=repmat(xbu_e',[1 N]);
 ybu2_e=repmat(ybu_e',[1 N]);
 xbv2_e=repmat(xbv_e',[1 N]);
 ybv2_e=repmat(ybv_e',[1 N]);
 
end
 
% HYCOM mask (determined by reading the first snapshot of SSH)
% in0: a list of masked points of the grid (XX,YY), presented as a 1d array
% iclosest = a 1d array of grid (XX,YY) indices of points that
% closest to the masked points from list in0
disp('determine hycom mask and interior cells closest to masked...');

if north_yes
 SSH=ncread(url_2d1,'ssh',[i1_n j1_n it1],[nx_n ny_n 1]);
 mask_hycom=ones(size(XX_n));
 mask_hycom(isnan(SSH))=0;
 [in0_n,iclosest_n]=find_closest(XX_n,YY_n,mask_hycom);
end

if south_yes
 SSH=ncread(url_2d1,'ssh',[i1_s j1_s it1],[nx_s ny_s 1]);
 mask_hycom=ones(size(XX_s));
 mask_hycom(isnan(SSH))=0;
 [in0_s,iclosest_s]=find_closest(XX_s,YY_s,mask_hycom);
end

if west_yes
 SSH=ncread(url_2d1,'ssh',[i1_w j1_w it1],[nx_w ny_w 1]);
 mask_hycom=ones(size(XX_w));
 mask_hycom(isnan(SSH))=0;
 [in0_w,iclosest_w]=find_closest(XX_w,YY_w,mask_hycom);
end

if east_yes
 SSH=ncread(url_2d1,'ssh',[i1_e j1_e it1],[nx_e ny_e 1]);
 mask_hycom=ones(size(XX_e));
 mask_hycom(isnan(SSH))=0;
 [in0_e,iclosest_e]=find_closest(XX_e,YY_e,mask_hycom);
end

for it=it_start:nt

 % Define the url directory
 if it==it_start
  kf=kf1;
  new_url=1;
  ihy=it1; % index of 1st ROMS rec. in the hycom time ser., found earlier
 end

 if new_url

  urldir=[urlhead urlends{kf}];
  url_2d=[urldir '2d'];
  url_uvel=[urldir 'uvel'];
  url_vvel=[urldir 'vvel'];
  url_temp=[urldir 'temp'];
  url_salt=[urldir 'salt'];  

  hycom_date=ncread(url_2d,'Date');
  nt_hy=length(hycom_date);

  new_url=0;

 end

 % read time:
 hycom_date_1=ncread(url_2d,'Date',[ihy],[1]);
 disp(' ');
 disp(int2str(hycom_date_1));
 day_h_1=datenum(int2str(hycom_date_1),'yyyymmdd');
 roms_day=day_h_1-datenum(date_ref)+day_ref; 

 %%% NORTH:

 if north_yes 
 
  %- read/interpolate:
  [SSH,U,V,Thycom,Shycom,ssh,u,v,T,S]=load_hycom_interp_2_roms...
     (url_2d, url_uvel, url_vvel, url_temp, url_salt,...
     XX_n, YY_n, z, i1_n, j1_n, in0_n, iclosest_n,ihy,...
     xbr_n, ybr_n, xbr2_n, ybr2_n, xbu2_n, ybu2_n, xbv2_n, ybv2_n,...
     zr_n, zu_n, zv_n);
 
  %- ubar,vbar
  ubar=get_ubar(u,Hz_u_n);
  vbar=get_ubar(v,Hz_v_n);
 
  %- Plot HYCOM and interpolated ROMS fields:
  cxs_u=[-0.5 0.5];
  cxs_v=[-0.5 0.5];
  cxs_T=[7 30];
  cxs_S=[36 37];
  xlims=[lon_rho(1,1) lon_rho(end,1)];
 
  kc=3;
  jo=floor(ny_n/2);
  Lon1=Lon(ii_n,jj_n(1)+jo-1);
  [zz,ll]=meshgrid(z,Lon1);
  
           
  %- output:
  ncwrite(fnameout,'ocean_time',roms_day,[it]);
  ncwrite(fnameout,'zeta_north',ssh,[1 it]);
  ncwrite(fnameout,'ubar_north',ubar,[1 it]);
  ncwrite(fnameout,'vbar_north',vbar,[1 it]);
  ncwrite(fnameout,'u_north',u,[1 1 it]);
  ncwrite(fnameout,'v_north',v,[1 1 it]);
  ncwrite(fnameout,'temp_north',T,[1 1 it]);
  ncwrite(fnameout,'salt_north',S,[1 1 it]);

 end
 
%%% SOUTH:

 if south_yes
 
  %- read/interpolate:
  [SSH,U,V,Thycom,Shycom,ssh,u,v,T,S]=load_hycom_interp_2_roms...
     (url_2d, url_uvel, url_vvel, url_temp, url_salt,...
     XX_s, YY_s, z, i1_s, j1_s, in0_s, iclosest_s,ihy,...
     xbr_s, ybr_s, xbr2_s, ybr2_s, xbu2_s, ybu2_s, xbv2_s, ybv2_s,...
     zr_s, zu_s, zv_s);
 
  %- ubar,vbar
  ubar=get_ubar(u,Hz_u_s);
  vbar=get_ubar(v,Hz_v_s);
 
  %- Plot HYCOM and interpolated ROMS fields:
  xlims=[lon_rho(1,1) lon_rho(end,1)];
 
  kc=1;
  jo=floor(ny_s/2);
  Lon1=Lon(ii_s,jj_s(1)+jo-1);
  [zz,ll]=meshgrid(z,Lon1);
 
  %- output:
  ncwrite(fnameout,'zeta_south',ssh,[1 it]);
  ncwrite(fnameout,'ubar_south',ubar,[1 it]);
  ncwrite(fnameout,'vbar_south',vbar,[1 it]);
  ncwrite(fnameout,'u_south',u,[1 1 it]);
  ncwrite(fnameout,'v_south',v,[1 1 it]);
  ncwrite(fnameout,'temp_south',T,[1 1 it]);
  ncwrite(fnameout,'salt_south',S,[1 1 it]);

 end 

%%% WEST:
 
 if west_yes
  
  %- read/interpolate:
  [SSH,U,V,Thycom,Shycom,ssh,u,v,T,S]=load_hycom_interp_2_roms...
     (url_2d, url_uvel, url_vvel, url_temp, url_salt,...
     XX_w, YY_w, z, i1_w, j1_w, in0_w, iclosest_w,ihy,...
     xbr_w, ybr_w, xbr2_w, ybr2_w, xbu2_w, ybu2_w, xbv2_w, ybv2_w,...
     zr_w, zu_w, zv_w);

  %- ubar,vbar
  ubar=get_ubar(u,Hz_u_w);
  vbar=get_ubar(v,Hz_v_w);
 
  %- Plot HYCOM and interpolated ROMS fields:
  xlims=[lat_rho(1,1) lat_rho(1,end)];
 
  kc=2;
  io=floor(nx_w/2);
  Lat1=Lat(ii_w(1)+io-1,jj_w);
  [zz,ll]=meshgrid(z,Lat1);
  
  %- output:
  ncwrite(fnameout,'zeta_west',ssh',[1 it]);
  ncwrite(fnameout,'ubar_west',ubar,[1 it]);
  ncwrite(fnameout,'vbar_west',vbar,[1 it]);
  ncwrite(fnameout,'u_west',u,[1 1 it]);
  ncwrite(fnameout,'v_west',v,[1 1 it]);
  ncwrite(fnameout,'temp_west',T,[1 1 it]);
  ncwrite(fnameout,'salt_west',S,[1 1 it]);

 end 

 %%% EAST:

 if east_yes
 
  %- read/interpolate:
  [SSH,U,V,Thycom,Shycom,ssh,u,v,T,S]=load_hycom_interp_2_roms...
     (url_2d, url_uvel, url_vvel, url_temp, url_salt,...
     XX_e, YY_e, z, i1_e, j1_e, in0_e, iclosest_e,ihy,...
     xbr_e, ybr_e, xbr2_e, ybr2_e, xbu2_e, ybu2_e, xbv2_e, ybv2_e,...
     zr_e, zu_e, zv_e);

  %- ubar,vbar
  ubar=get_ubar(u,Hz_u_e);
  vbar=get_ubar(v,Hz_v_e);
 
  %- Plot HYCOM and interpolated ROMS fields:
  xlims=[lat_rho(1,1) lat_rho(1,end)];
 
  kc=4;
  io=floor(nx_e/2);
  Lat1=Lat(ii_e(1)+io-1,jj_e);
  [zz,ll]=meshgrid(z,Lat1);

  %- output:
  ncwrite(fnameout,'zeta_east',ssh',[1 it]);
  ncwrite(fnameout,'ubar_east',ubar,[1 it]);
  ncwrite(fnameout,'vbar_east',vbar,[1 it]);
  ncwrite(fnameout,'u_east',u,[1 1 it]);
  ncwrite(fnameout,'v_east',v,[1 1 it]);
  ncwrite(fnameout,'temp_east',T,[1 1 it]);
  ncwrite(fnameout,'salt_east',S,[1 1 it]);

 end % east_yes
 
 ihy=ihy+1;

 if ihy>nt_hy
  new_url=1;
  kf=kf+1;
  ihy=1;
 end  
 
end  % for it=1:nt

