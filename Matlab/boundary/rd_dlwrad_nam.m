clear all;
wrkdir='/home/server/student/homes/ipasmans/PhD/projects/boundary/';
ldir=[wrkdir 'matlib/'];
path(ldir,path);
path([ldir 'AIR_SEA'],path);
path([ldir 'FUNCTIONS'],path);
eval(['load ' '/home/aruba/vol1/serofeev/OW/ETA/limits.dat']);
lat1=limits(3);lat2=limits(4);lon1=limits(1);lon2=limits(2);
%
date_ref='01-Jan-2005'; % ref date and time in roms application
day_ref=0;
timeunits=['days since ' datestr(date_ref,'yyyy-mm-dd HH:MM:SS')];
%
% dayn=floor(now); %REMOVE -1!
% today=datestr(dayn);
% yester=datestr(dayn-1); % data available until yesterday
% date0=datestr(dayn-8);
roms_dates=datestr({'2011-04-01','2011-04-02'});
%
Gname=[wrkdir 'Prm/grd_ow2km_roms.nc'];
mask_rho=ncread(Gname,'mask_rho');
lon=ncread(Gname,'lon_rho');
lat=ncread(Gname,'lat_rho');
%
[Mr,Lr]=size(mask_rho); % same as [eta_rho,xi_rho]
eta_rho=Lr;xi_rho=Mr;

% Example URL:
% url='http://nomads.ncdc.noaa.gov:80/dods/NCEP_NAM/201304/20130401/nam_218_20130401_0000_fff'
urlhead='http://nomads.ncdc.noaa.gov:80/dods/NCEP_NAM/';


dstr1=datestr(roms_dates{1},'yyyymmdd');
dstr2=datestr(roms_dates{1},'yyyymmdd');
day1=datenum(roms_dates{1});
day2=datenum(roms_dates{2});

days=[day1:day2];
ndays=length(days);

% read NAM lon, lat, choose NAM indices to read from:
url=[urlhead dstr1(1:6) '/' dstr1 '/nam_218_' dstr1 '_0000_fff'];

     
nam_lon=ncread(url,'lon');
nam_lat=ncread(url,'lat');
ii=findin(nam_lon,[lon1 lon2]+[-0.2,0.2]);
jj=findin(nam_lat,[lat1 lat2]+[-0.2,0.2]);
i1=ii(1);
j1=jj(1);
nx=length(ii);
ny=length(jj);


[dllat,dllon]=meshgrid(nam_lat(jj),nam_lon(ii));

% Read the list of unavailable dates, make the list of avail. times 
% and the rec length 
tdays_avail=[day1:6/24:day2]';
ndays_avail=length(tdays_avail);
a=load('nam_available_url.txt');
na=size(a,1); % a(:,1:3) - yyyy mm dd, a(:,4:7) - 0 6 12 18, if avail. 
b=a(:,4:7);   % nan if unavail.
list_days=datenum(a(:,1),a(:,2),a(:,3));
hr=[0 6 12 18];
list_days=repmat(list_days,[1 4])+repmat(hr,[na 1])/24;
list_days=reshape(list_days',[na*4 1]);
inans=find(isnan(reshape(b',[na*4 1])));
nan_days=list_days(inans);
in=findin(nan_days,[day1 day2]);
nan_days=nan_days(in);
if ~isempty(nan_days)
 for ka=1:length(nan_days)
  tdays_avail(find(tdays_avail==nan_days(ka)))=NaN;   
 end
end
tdays_avail=tdays_avail(~isnan(tdays_avail));
nrec=diff([tdays_avail;day2])/0.125; % the length of the record to read
nrec(end)=nrec(end)+1;               % from each avail. url                        

tdays=[day1:3/24:day2];  % time ser, proposal
[nl,ml]=size(dllat);
dlwrad_ts=zeros(nl,ml,1000);
time_dl=[];irec=1;Td=[];
for k=1:ndays_avail
 dstr=datestr(tdays_avail(k),'yyyymmdd HH');
 hr=[dstr(end-1:end) '00'];
 dstr=dstr(1:8);
 %url=[urlhead dstr(1:6) '/' dstr '/nam_218_' dstr '_hh00_000'];
 
 url=[urlhead dstr(1:6) '/' dstr '/nam_218_' dstr '_' hr '_fff'];
 
 
 nr=nrec(k);
 for it=1:nr
     
  if k==1 && it==1
   land=ncread(url,'landsfc',[i1 j1 1],[nx ny 1]);
   land(land>0)=1;
   nam_mask=1-land;
   [in0,iclosest]=find_closest(dllon,dllat,nam_mask);      
      http://nomads.ncdc.noaa.gov/data/namanl/201104/20110401/namanl_218_20110401_0000_000.grb
  end
  
  t=ncread(url,'time',[it],[1])+datenum(1,1,1)-2;
  disp(datestr(t));
  
  url=['http://nomads.ncdc.noaa.gov/data/namanl/',datestr(t,'yyyymm'),...
      '/',datestr(t,'yyyyymmdd'),'/namanl_218_',datestr(t,'yyyymmdd'),'_',...
      
      
      http://nomads.ncdc.noaa.gov/data/namanl/201104/20110401/namanl_218_20110401_0000_000.grb
  
  dlwrad=ncread(url,'dlwrfsfc',[i1 j1 it],[nx ny 1]);  
  dlwrad_ts(:,:,irec)=dlwrad;
  troms=t-datenum(date_ref)+day_ref;
  time_dl=[time_dl;troms];Td=[Td;t];  
  irec=irec+1;
 end 
end
nt=irec-1;
dlwrad_ts=dlwrad_ts(:,:,1:nt);
save DLW.mat dlwrad_ts time_dl Td dllon dllat;
break
figure(1);
ndays_avail=length(tdays_avail);
for irec=1:nt
  clf;
  pcolor(dllon,dllat,dlwrad_ts(:,:,irec));shading flat;
  caxis([0 500]);
  title(datestr(Td(irec)));  
  waitforbuttonpress;
end
