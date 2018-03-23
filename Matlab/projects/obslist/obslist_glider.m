function [ ] = obslist_glider(grdfile,dataDir,outfile,tLim,mask)
%obslist_glider Create obslist based on glider data
%   oblist_glider(grdfile,dataDir,outfile,tLim)

addpath('/home/server/student/homes/ipasmans/Matlab/roms/seawater/'); 

flag_temp=true; 
flag_salt=true; 
flag_thin=true; 

%% observational uncertainty

R0=1025; Scoef=7.5e-4; Tcoef=1.7e-4; drho=1e-2; 
sig_T=load('/home/server/student/homes/ipasmans/Matlab/projects/obslist/interp_error_temp_h1000v2.mat'); 
sig_T.fun=@(x) .5*interp1([-100e3,sig_T.z(1:end-1),100e3],sqrt([sig_T.Var(1),sig_T.Var',sig_T.Var(end)]),x); 
sig_T.min=drho/R0/Tcoef; 
sig_S=load('/home/server/student/homes/ipasmans/Matlab/projects/obslist/interp_error_salt_h1000v2.mat'); 
sig_S.fun=@(x) .5*interp1([-100e3,sig_S.z(1:end-1),100e3],sqrt([sig_S.Var(1),sig_S.Var',sig_S.Var(end)]),x); 
sig_S.min=drho/R0/Scoef; 


%% read grid

grd.lon=ncread(grdfile,'lon_rho'); 
grd.lat=ncread(grdfile,'lat_rho'); 
grd.mask=ncread(grdfile,'mask_rho');
grd.zw=ncread(grdfile,'z0_w'); 
grd.h=ncread(grdfile,'h'); 

if isempty(mask)
    mask=grd.mask; 
end

%% read glider data

con=dir(dataDir); 
obs=struct('t',[],'temp',[],'salt',[],'sig_temp',[],'sig_salt',[],'depth',[],'lon',[],'lat',[]); 

for iFile=1:length(con)
   if ~strncmpi(con(iFile).name,'NH',2); continue; end
   
   %glider file
   fileIn=fullfile(dataDir,con(iFile).name); 
   display(fileIn); 
   data=load(fileIn); 
   data=data.gridded_data; 

   %read time
   obs1.t=data.Timeinsec/24/3600+datenum('1970-01-01'); 
   if ~any(obs1.t>=tLim(1) & obs1.t<=tLim(2)); continue; end
   s=size(data.Temp); 
   obs1.t=repmat(obs1.t,[1,s(2)]); 
   
   %read lat/lon
   obs1.lon=mod(repmat(data.Lon,[1,s(2)]),360); 
   obs1.lon(obs1.lon>180)=obs1.lon(obs1.lon>180)-360;
   obs1.lat=repmat(data.Lat,[1,s(2)]); 
   
   %depth
   obs1.depth=sw_dpth(data.Press,obs1.lat); 
   
   %temperature 
   obs1.temp=data.Temp;
   obs1.sig_temp=data.Temp_std; 
   
   %salinity
   obs1.salt=data.Salt; 
   obs1.sig_salt=data.Salt_std; 
   
   %collect data
   inmask=interp2(grd.lon',grd.lat',mask',obs1.lon(:),obs1.lat(:)); 
   inan=obs1.t(:)>=tLim(1)&obs1.t(:)<=tLim(2)&~isnan(obs1.depth(:))&~isnan(obs1.temp(:))&~isnan(obs1.salt(:))...
       &obs1.depth(:)>0&inmask==1; 
   
   obs.t=[obs.t;obs1.t(inan)];
   obs.temp=[obs.temp;obs1.temp(inan)]; 
   obs.sig_temp=[obs.sig_temp;obs1.sig_temp(inan)];
   obs.salt=[obs.salt;obs1.salt(inan)]; 
   obs.sig_salt=[obs.sig_salt;obs1.sig_salt(inan)]; 
   obs.depth=[obs.depth;obs1.depth(inan)]; 
   obs.lon=[obs.lon;obs1.lon(inan)]; 
   obs.lat=[obs.lat;obs1.lat(inan)]; 
   
   %add other uncertainties
   %Fore exp35 and later disable sig_T.fun and sig_S.fun
   obs.sig_temp=sqrt( obs.sig_temp.^2 ...
       +sig_T.min^2*ones(size(obs.sig_temp))...
       +0*sig_T.fun(-obs.depth).^2 ); 
   obs.sig_salt=sqrt( obs.sig_salt.^2 ...
       +sig_S.min^2*ones(size(obs.sig_temp))...
       +0*sig_S.fun(-obs.depth).^2 ); 
   
end

%% Thin

%Thin observations such that there is ~1 per layer
if flag_thin
    display('Thin glider'); 
    
    obs1=obs; clear obs;
    obs=struct('t',[],'temp',[],'salt',[],'sig_temp',[],'sig_salt',[],'depth',[],'lon',[],'lat',[]); 
    [ll,ll2uni,uni2ll]=unique([obs1.lon,obs1.lat],'rows');
    
    
    obs1.flag=zeros(size(obs1.lon)); 
    for i1=1:size(ll,1)
       %Location in grid indices
       [~,ilon]=min(abs(grd.lon(:,1)-ll(i1,1))); 
       [~,ilat]=min(abs(grd.lat(1,:)-ll(i1,2))); 
       
      
       for i3=1:size(grd.zw,3)-1
           in=obs1.flag==0 & uni2ll==i1 & ...
               -obs1.depth>=grd.zw(ilon,ilat,i3) & -obs1.depth<=grd.zw(ilon,ilat,i3+1);
           if ~any(in); continue; end
           
           obs.t=[obs.t;mean(obs1.t(in))]; 
           obs.lon=[obs.lon;mean(obs1.lon(in))]; 
           obs.lat=[obs.lat;mean(obs1.lat(in))]; 
           obs.depth=[obs.depth;mean(obs1.depth(in))]; 
           obs.temp=[obs.temp;mean(obs1.temp(in))]; 
           obs.salt=[obs.salt;mean(obs1.salt(in))]; 
           obs.sig_temp=[obs.sig_temp;hypot(std(obs1.temp(in),1),rms(obs1.sig_temp(in)))]; 
           obs.sig_salt=[obs.sig_salt;hypot(std(obs1.salt(in),1),rms(obs1.sig_salt(in)))]; 
           obs1.flag(in)=1; 
       end
    end
    
    
end

%% create outfile if necessary

%create output
if exist(outfile,'file')
 Ktotal=ncdim1(outfile,'K');
else
 % create the output file (obslist):
 nccreate(outfile,'obs','Dimensions',{'K' Inf},'format','classic');
 nccreate(outfile,'sig_d','Dimensions',{'K' Inf});
 nccreate(outfile,'type','Dimensions',{'K' Inf},'Datatype','int32');
 nccreate(outfile,'lon','Dimensions',{'K' Inf});
 nccreate(outfile,'lat','Dimensions',{'K' Inf});
 nccreate(outfile,'z','Dimensions',{'K' Inf});
 nccreate(outfile,'time','Dimensions',{'K' Inf});

 Ktotal=0;
end

%% write data to file
 

display('write temperature data to file');
K=length(obs.t);
if K>0

    %temp
    if flag_temp
    ncwrite(outfile,'obs',obs.temp,Ktotal+1);
    ncwrite(outfile,'sig_d',obs.sig_temp,Ktotal+1);
    ncwrite(outfile,'type',6*ones(K,1),Ktotal+1);
    ncwrite(outfile,'lon',obs.lon,Ktotal+1);
    ncwrite(outfile,'lat',obs.lat,Ktotal+1);
    ncwrite(outfile,'z',obs.depth,Ktotal+1);
    ncwrite(outfile,'time',(obs.t-tLim(1))*24*3600,Ktotal+1);
    Ktotal=Ktotal+K; 
    end
    
    %salt
    if flag_salt
    ncwrite(outfile,'obs',obs.salt,Ktotal+1);
    ncwrite(outfile,'sig_d',obs.sig_salt,Ktotal+1);
    ncwrite(outfile,'type',7*ones(K,1),Ktotal+1);
    ncwrite(outfile,'lon',obs.lon,Ktotal+1);
    ncwrite(outfile,'lat',obs.lat,Ktotal+1);
    ncwrite(outfile,'z',obs.depth,Ktotal+1);
    ncwrite(outfile,'time',(obs.t-tLim(1))*24*3600,Ktotal+1);
    Ktotal=Ktotal+K;
    end
end

display(sprintf('%d observations written',2*K));
end

