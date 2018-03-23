function obs=obslist_l3u(goesfile,grdfile,outfile,tLim,mask)

display(sprintf('%s: Starting obslist_l3u',datestr(now(),'dd mmm HH:MM'))); 

flag_thin=true; 

%Based std SST adjacent grid points for Exp26f_2353 and |dSST|<.45:
% dSST/dxi:0.0674 C; dSST/deta: 0.0551 C; dSST/dt: 0.053 C based on std
% dSST/dxi:0.0639 C; dSST/deta: 0.0499 C; dSST/dt: 0.0497 C for xi 1-240
% dSST/dxi:0.0251 C; dSST/deta: 0.0232 C; dSST/dt: .0232 C based on Gaussian fit 
sig0=0.0; 

% create the output file (obslist):
if ~exist(outfile,'file')
    Ktotal=0;
    nccreate(outfile,'obs','Dimensions',{'K' Inf},'format','classic');
    nccreate(outfile,'sig_d','Dimensions',{'K' Inf});
    nccreate(outfile,'type','Dimensions',{'K' Inf},'Datatype','int32');
    nccreate(outfile,'lon','Dimensions',{'K' Inf});
    nccreate(outfile,'lat','Dimensions',{'K' Inf});
    nccreate(outfile,'z','Dimensions',{'K' Inf});
    nccreate(outfile,'time','Dimensions',{'K' Inf});
else
    Ktotal=length(ncread(outfile,'obs')); 
end

%read grid
grd.lon=ncread(grdfile,'lon_rho'); 
grd.lat=ncread(grdfile,'lat_rho'); 
grd.mask=ncread(grdfile,'mask_rho'); 
grd.lon2=grd.lon(1:2:end,1:2:end); 
grd.lat2=grd.lat(1:2:end,1:2:end); 


%% read data

obs=struct('lon',[],'lat',[],'t',[],'temp',[],'val',[],'sig',[],'errorMean',[]); 
for iFile=1:length(goesfile)
    %load
    data=load(goesfile{iFile}); data=data.obs; 
   
    for it=1:length(data)
       if data(it).t>=tLim(1) && data(it).t<=tLim(2) && ~ismember(data(it).t,obs.t)
           
           data(it).lat=reshape(data(it).lat,[],1); 
           data(it).lon=reshape(data(it).lon,[],1); 
           
          in= ...
              interp2(grd.lon',grd.lat',mask',data(it).lon,data(it).lat)>=1; 
          in=in & data(it).quality==5 & ~isnan(data(it).std) & ~isnan(data(it).bias); 
          obs.lon=[obs.lon;data(it).lon(in)]; 
          obs.lat=[obs.lat;data(it).lat(in)]; 
          obs.t=[obs.t;data(it).t*ones(sum(in),1)]; 
          obs.temp=[obs.temp;data(it).sst(in)-1*data(it).bias(in)]; 
          obs.errorMean=[obs.errorMean;data(it).bias(in)]; 
          obs.sig=[obs.sig;data(it).std(in)]; 
       end
    end
    
end

%% Thin data

if flag_thin
    dt=3/24; %time in days
    
    %Bin data in space
    display('Binning in space'); 
    obs1=obs; 
    for i1=1:length(obs.lon)
        obs1.ilon(i1,1)=find(obs.lon(i1)>=grd.lon2(1:end-1,1)&obs.lon(i1)<=grd.lon2(2:end,1),1); 
        obs1.ilat(i1,1)=find(obs.lat(i1)>=grd.lat2(1,1:end-1)&obs.lat(i1)<=grd.lat2(1,2:end),1); 
    end
    
    clear obs; 
    [llUni,ll2uni,uni2ll]=unique([obs1.ilon,obs1.ilat,obs1.t],'rows'); 
    for i1=1:size(llUni,1)
       in=uni2ll==i1; 
       obs.lon(i1,1)=mean(obs1.lon(in)); 
       obs.lat(i1,1)=mean(obs1.lat(in)); 
       obs.t(i1,1)=mean(obs1.t(in)); 
       obs.temp(i1,1)=mean(obs1.temp(in)); 
       obs.errorMean(i1,1)=mean(obs1.errorMean(in)); 
       obs.sig(i1,1)=hypot(std(obs1.temp(in),1),rms(obs1.sig(in))); 
       obs.ilon(i1,1)=mean(obs1.ilon(in)); 
       obs.ilat(i1,1)=mean(obs1.ilat(in)); 
    end
    display(sprintf('Thinning factor: %.2f',length(obs.lon)/length(obs1.lon))); 
    
    %Find times on which SST observations are available
    display('Thinning in time'); 
    tUni=unique(obs.t); 
    tUni=sort(tUni,'descend'); 
    [llUni,ll2uni,uni2ll]=unique([obs.ilon,obs.ilat],'rows'); 
    
    obs.last=Inf*ones(size(obs.lon)); 
    obs.flag=zeros(size(obs.lon)); 
    for it=1:length(tUni)
        in= obs.t==tUni(it) & obs.last>=tUni(it)+dt; 
        obs.flag(in)=1; 
        
        in=ismember(uni2ll,uni2ll(in)); 
        obs.last(in)=tUni(it); 
    end  
    display(sprintf('Thinning factor: %.2f',sum(obs.flag==1)/length(obs.lon))); 
    
    in=obs.flag==1; 
    obs.lon=obs.lon(in); 
    obs.lat=obs.lat(in); 
    obs.t=obs.t(in); 
    obs.temp=obs.temp(in); 
    obs.errorMean=obs.errorMean(in); 
    obs.sig=obs.sig(in); 
    obs.sig=hypot(obs.sig,obs.errorMean); 
    obs.ilon=obs.ilon(in); 
    obs.ilat=obs.ilat(in); 
    

end

%% make ready for output

%obs.sig=sqrt(obs.sig.^2+0*nanmax(abs(obs.errorMean))+sig0^2); 
obs.sig= sqrt(sig0^2+obs.sig.^2+0*obs.errorMean.^2); 
obs.type=4*ones(size(obs.t)); 
obs.z=zeros(size(obs.t)); 
obs.t=(obs.t-tLim(1))*24*3600; 

%remove NaN
in=~isnan(obs.temp) ; 
obs.lon=obs.lon(in); 
obs.lat=obs.lat(in); 
obs.type=obs.type(in); 
obs.t=obs.t(in); 
obs.temp=obs.temp(in); 
obs.sig=obs.sig(in); 
obs.z=obs.z(in); 

%write to output
info=ncinfo(outfile,'type'); 
Ktotal=info.Size(1); 
K=length(obs.t); 
if K>0
    ncwrite(outfile,'obs',obs.temp,Ktotal+1);
    ncwrite(outfile,'sig_d',obs.sig,Ktotal+1);
    ncwrite(outfile,'type',obs.type,Ktotal+1);
    ncwrite(outfile,'lon',obs.lon,Ktotal+1);
    ncwrite(outfile,'lat',obs.lat,Ktotal+1);
    ncwrite(outfile,'z',obs.z,Ktotal+1);
    ncwrite(outfile,'time',obs.t,Ktotal+1);
end


end


