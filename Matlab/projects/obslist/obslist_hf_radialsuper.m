function obs1=obslist_hf_radialsuper(grdfile,obsfile,tide,outfile,tLim,mask)
% obslist_hf(grdfile,obs_dir,outfile,tLim,tOut,mask)
display(sprintf('%s: Starting obslist_hf_radialsuper',datestr(now(),'dd mmm HH:MM'))); 
addpath('/home/server/student/homes/ipasmans/Matlab/roms/t_tide'); 
addpath('/home/server/student/homes/ipasmans/Matlab/figures/'); 

flag_tide=false; 
flag_filter=false;  


%Create output
typeo=[]; sig_d=[]; val=[];
lon=[]; lat=[]; time=[]; theta=[]; 


%Std daily-average daily signal 2011-4-1 0200-2011-9-29 2300: u:.0018/v:.0022 m/s
%Std daily-averaged HF signal: u:0.0149/v:.0146 m/s
%Std daily-averaged HF signal minus tide: u:0.0148/v:.0144 m/s
%For grid xi 1-240 eta 1-Inf
%Std daily-averaged du_daily/dxi on Exp26f_2353: 0.0111/v
%Std daily-avraged du_daily/deta on Exp26f_2353: 0.0166 m/s
%Std daily-averaged dv_daily/dxi on Exp26f_2353: 0.0210/v
%Std daily-avraged dv_daily/deta on Exp26f_2353: 0.0112 m/s
sig0=0;

%% Read grid

lon_r=ncread(grdfile,'lon_rho'); 
lat_r=ncread(grdfile,'lat_rho'); 
lon_u=ncread(grdfile,'lon_u');
lat_u=ncread(grdfile,'lat_u');
lon_v=ncread(grdfile,'lon_v');
lat_v=ncread(grdfile,'lat_v');
mask_r=ncread(grdfile,'mask_rho'); 
mask_u=ncread(grdfile,'mask_u');
mask_v=ncread(grdfile,'mask_v');
deg2m=[cos(deg2rad(mean(lat_r(:)))),1]*2*pi*earthRadius/360;


%% Read and process files

obs1=struct('lon',[],'lat',[],'z',[],'t',[],'val',[],'sig',[],'dir',[],'type',[]);

%read files
for iFile=1:length(obsfile)
    fname=obsfile{iFile};
    
    %Load data
    obs=load(fname); obs=obs.obs;
    
    %Select suitable measurements
    mu=interp2(lon_u',lat_u',mask_u',obs.lon,obs.lat);
    mv=interp2(lon_v',lat_v',mask_v',obs.lon,obs.lat);
    mr=interp2(lon_r',lat_r',mask',obs.lon,obs.lat);
    
    in=find(mu==1 & mv==1 & mr==1 );
    if ~any(in); continue; end
    obs.lon=obs.lon(in); obs.lat=obs.lat(in);
    obs.vr=obs.vr(in,:); obs.theta=obs.theta(in);
    obs.sig_space=obs.sig_space(in,:);
    obs.sig_time=obs.sig_time(in,:);
    
    %Radial locations
    obs.x=deg2rad(obs.lon-obs.lon0)*cos(deg2rad(obs.lat0))*earthRadius;
    obs.y=deg2rad(obs.lat-obs.lat0)*earthRadius;
    obs.r=hypot(obs.x,obs.y);
    obs.phi=mod(rad2deg(atan2(obs.x,obs.y)),360);
    
    %Replace spatial error with location average
    inCount=0;
    for i1=1:size(obs.sig_space,1)
        in=obs.t>=tLim(1) & obs.t<=tLim(2) & ~isnan(obs.sig_space(i1,:));
        if ~any(in)
            obs.sig_space(i1,isnan(obs.sig_space(i1,:)))=nanmedian(obs.sig_space(i1,:));
        else
            inCount=inCount+1;
            obs.sig(i1,isnan(obs.sig_space(i1,:)))=nanmedian(obs.sig_space(i1,in));
        end
    end
    display(sprintf('Space error fraction %s: %.2f',fname,inCount/size(obs.sig_space,1)));
    
    %% Filter
    
    %LP-filter
    if flag_filter
        in=obs.t>=tLim(1)-diff(tLim) & obs.t<=tLim(2)+diff(tLim);
        
        %Select times in extended window
        obs.t=obs.t(in);
        obs.vr=obs.vr(:,in);
        obs.sig_time=obs.sig_time(:,in);
        obs.sig_space=obs.sig_space(:,in);
        
        %Fill in missing observations
        tf=[min(obs.t):1/24:max(obs.t)];
        obs.vr=fill_interp(obs.t,obs.vr,tf,3/24,2);
        obs.vr=filter_lanczos(obs.vr,1,1/40,60,2);
        obs.sig_time=fill_interp(obs.t,obs.sig_time,tf,3/24,2);
        obs.sig_space=fill_interp(obs.t,obs.sig_space,tf,3/24,2);
        obs.t=tf;
    else
        in=obs.t>=tLim(1) & obs.t<=tLim(2);
        
        %Select times in window
        obs.t=obs.t(in);
        obs.vr=obs.vr(:,in);
        obs.sig_time=obs.sig_time(:,in);
        obs.sig_space=obs.sig_space(:,in);
    end
    
    if flag_tide
        %Calculate tidal components
        tidecon=zeros(size(tide.freq,1),8,length(obs.lon));
        for i3=1:length(tide.freq)
            tidecon(i3,1,:)=tide.major_interp{i3}(obs.lon,obs.lat);
            tidecon(i3,3,:)=tide.minor_interp{i3}(obs.lon,obs.lat);
            tidecon(i3,5,:)=tide.dir_interp{i3}(obs.lon,obs.lat);
            tidecon(i3,7,:)=tide.phase_interp{i3}(obs.lon,obs.lat);
            tidecon(i3,5:8,:)=rad2deg(angle(tidecon(i3,5:8,:)));
        end
        obs.tidecon=tidecon;
        
        %Calculate tide
        for i1=1:length(obs.lon)
            if any(isnan(tidecon(:,:,i1))); continue; end
            tt=struct('name',tide.name,'freq',tide.freq,'tidecon',squeeze(tidecon(:,:,i1)),'type','nodal');
            vtide=t_predic(obs.t,tt,'latitude',obs.lat(i1));
            obs.vtide(i1,:)=real(vtide)*sin(deg2rad(obs.theta(i1)))+imag(vtide)*cos(deg2rad(obs.theta(i1)));
        end
        
        obs.vr=obs.vr+obs.vtide;
    else
        obs.vr=obs.vr;
    end
    
    %% Cluster observations
    
    phi0=min(obs.phi);
    obs.phi=mod(obs.phi-phi0,360); 
    
    
    %Cluster in time
    dt=1; tDay=[tLim(1)+.5*dt:dt:tLim(2)];
    for it=1:length(tDay)
        in2=obs.t>=tDay(it)-.5*dt & obs.t<=tDay(it)+.5*dt;
        if sum(in2)<=23; continue; end
        
        %Cluster in location
        dr=5e3; r=[1.5*dr:dr:max(obs.r)];
        for ir=1:length(r) 
            dphi=360/floor(2*pi*r(ir)/dr);
            phi=[.5*dphi:dphi:360];
            for iphi=1:length(phi)
                %Velocities within time and location frame
                in1=obs.phi>=phi(iphi)-.5*dphi & obs.phi<phi(iphi)+.5*dphi & ...
                    obs.r>=r(ir)-.5*dr & obs.r<r(ir)+.5*dr;
                if ~any(in1); continue; end
                
                %Force radial directions in bin to be similar
                if max(mod(obs.theta(in1)-min(obs.theta(in1)),180))>min(dphi,9) || ...
                        max(mod(max(obs.theta(in1))-obs.theta(in1),180))>min(dphi,9)
                    continue; 
                end
                
                %Velocities in Cartesian coordinates
                u1=bsxfun(@times,obs.vr(in1,in2),sin(deg2rad(obs.theta(in1))));
                v1=bsxfun(@times,obs.vr(in1,in2),cos(deg2rad(obs.theta(in1))));
                nnan=sum(~isnan(u1),1); 
                u1=nanmean(u1,1); v1=nanmean(v1,1);
                
                %Uncertainties in measurement
                if sum(nnan==0)>2; continue; end
                u1=fill_nan(u1); v1=fill_nan(v1); 
                sigspace1=nansum(obs.sig_space(in1,in2).^2,1)./nnan.^2;
                sigtime1=nansum(obs.sig_time(in1,in2).^2,1)./nnan.^2;
                
                
                %Time average and space average
                nnan=sum(nnan~=0); 
                sigtime1=sqrt( nansum(sigtime1)/nnan^2 ); 
                sigspace1=sqrt( nansum(sigspace1)/nnan );
                u1=trapz(obs.t(in2),u1)/(max(obs.t(in2))-min(obs.t(in2))); 
                v1=trapz(obs.t(in2),v1)/(max(obs.t(in2))-min(obs.t(in2))); 
                
                %To output
                obs1.lon=[obs1.lon;mean(obs.lon(in1))];
                obs1.lat=[obs1.lat;mean(obs.lat(in1))];
                obs1.z=[obs1.z;0];
                obs1.t=[obs1.t;tDay(it)];
                obs1.val=[obs1.val;hypot(u1,v1)];
                obs1.sig=[obs1.sig;norm([sig0,sigtime1,sigspace1])];
                obs1.dir=[obs1.dir;atan2(u1,v1)];
                obs1.type=[obs1.type;8]; 
            end
        end
    end
    
end

%Adjust direction
in=cos(obs1.dir)<0; 
obs1.val(in)=-obs1.val(in); 
obs1.dir(in)=obs1.dir(in)+pi; 
obs1.dir=mod(rad2deg(obs1.dir),360); 

%Adjust time
obs1.t=(obs1.t-tLim(1))*24*3600; 

%Sort by time and remove NaNs
val1=sortrows([obs1.t,obs1.lon,obs1.lat,obs1.z,obs1.val,obs1.sig,obs1.dir,obs1.type]); 
val1=val1(~any(isnan(val1),2),:); 

obs1.t=val1(:,1); obs1.lon=val1(:,2); obs1.lat=val1(:,3); obs1.z=val1(:,4); 
obs1.val=val1(:,5); obs1.sig=val1(:,6); obs1.dir=val1(:,7); obs1.type=val1(:,8); 
 
%write to output
Ktotal=length(ncread(outfile,'obs')); 
if ~isempty(obs1.val)>0
    ncwrite(outfile,'obs',obs1.val,Ktotal+1);
    ncwrite(outfile,'sig_d',obs1.sig,Ktotal+1);
    ncwrite(outfile,'type',obs1.type,Ktotal+1);
    ncwrite(outfile,'lon',obs1.lon,Ktotal+1);
    ncwrite(outfile,'lat',obs1.lat,Ktotal+1);
    ncwrite(outfile,'z',obs1.z,Ktotal+1);
    ncwrite(outfile,'time',obs1.t,Ktotal+1);
    ncwrite(outfile,'dir',obs1.dir,Ktotal+1);
end


end

function val=fill_nan(val)
    %FILL_NAN fill in NaN values
    
    in=~isnan(val); 
    t=[1:length(val)]; 
    val=interp1(t(in),val(in),t); 
    j1=find(in,1,'first'); j2=find(in,1,'last'); 
    val(1:j1)=val(j1); val(j2:end)=val(j2); 
end
