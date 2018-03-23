%% Create obslist files
clc; clear all; 

%time 
tWin=3; %length window in days
dateref=datenum('2005-01-01'); %reference date model
tOut=[2380:3:2380]+dateref; %2423+dateref; 
%tOut=[2329:3:2461]+dateref; 
%grid file
grdfile='/home/aruba/vol2/ipasmans/Exp/Prm/grd_ow2km_r13lin_mix.nc'; 
%directory to which output will be written
outdir='/home/aruba/vol2/ipasmans/Exp/Obs/Exp37_2380'; 
%directory that contains observations
obsDir='/home/server/student/homes/ipasmans/Obs/';
%If sig_bal temp level falls below mask_level no observations are
%assimilated at the location
mask_level=1.15*.9/1.2; 

%% Create mask

if ~exist(outdir,'dir'); mkdir(outdir); end

sigfile='/home/aruba/vol2/ipasmans/Exp/Prm/per_sig_bal.nc'; 
mask=ncread(sigfile,'temp',[1,1,40,1],[Inf,Inf,1,1]);  mask=squeeze(mask); 
mask=mask>=mask_level;
mask=double(mask); 

grd.lon=ncread(grdfile,'lon_rho'); 
grd.lat=ncread(grdfile,'lat_rho'); 

%Remove points close to wall
in=abs(grd.lon-min(grd.lon(:))).*cos(deg2rad(grd.lat))*2*pi*earthRadius/360<50e3; 
mask(in)=0; 
in=abs(grd.lat-max(grd.lat(:)))*2*pi*earthRadius/360<50e3; 
mask(in)=0; 
in=abs(grd.lat-min(grd.lat(:)))*2*pi*earthRadius/360<50e3; 
mask(in)=0; 

%Remove points close to CR
in=abs(grd.lon+124.0407).*cos(deg2rad(grd.lat))*2*pi*earthRadius/360<10e3 & ...
    abs(grd.lat-46.2297)*2*pi*earthRadius/360<10e3;
mask(in)=0; 

%% Read tides

display('Generating interpolants for tide'); 

%Horizontal tide
tideUV=load('/home/server/student/homes/ipasmans/Matlab/projects/tide_analysis/tideUV40_Exp26fstride1.mat');
tideUV.phase=exp(1i*deg2rad(tideUV.phase));
tideUV.dir=exp(1i*deg2rad(tideUV.dir));

%Interpolants
in=~isnan(tideUV.major(:,:,1));
in=reshape(in,[],1);
for i3=1:size(tideUV.major,3)
    val1=squeeze(tideUV.major(:,:,i3));
    tideUV.major_interp{i3}=scatteredInterpolant(tideUV.lon(in),tideUV.lat(in),val1(in));
    val1=squeeze(tideUV.minor(:,:,i3));
    tideUV.minor_interp{i3}=scatteredInterpolant(tideUV.lon(in),tideUV.lat(in),val1(in));
    val1=squeeze(tideUV.dir(:,:,i3));
    tideUV.dir_interp{i3}=scatteredInterpolant(tideUV.lon(in),tideUV.lat(in),val1(in));
    val1=squeeze(tideUV.phase(:,:,i3));
    tideUV.phase_interp{i3}=scatteredInterpolant(tideUV.lon(in),tideUV.lat(in),val1(in));
end

%Vertical tide
tideZ=load('/home/server/student/homes/ipasmans/Matlab/projects/tide_analysis/tide_Exp26fstride1.mat');
in=~isnan(tideZ.val(:,:,1)); 
in=reshape(in,[],1); 
for i3=1:size(tideZ.val,3)
   val1=squeeze(tideZ.val(:,:,i3)); 
   tideZ.major_interp{i3}=scatteredInterpolant(tideZ.lon(in),tideZ.lat(in),abs(val1(in))); 
   tideZ.phase_interp{i3}=scatteredInterpolant(tideZ.lon(in),tideZ.lat(in),val1(in)./abs(val1(in))); 
end

%% generate obs file for each time window

for it=1:length(tOut)
    display(sprintf('Making obslist for %s',datestr(tOut(it)))); 
    
    %Name output
    obsfile=fullfile(outdir,sprintf('obslist%4d.nc',tOut(it)-dateref));
    if exist(obsfile,'file')
        %continue;
        delete(obsfile);
    end
    nccreate(obsfile,'obs','Dimensions',{'K' Inf},'format','classic');
    nccreate(obsfile,'sig_d','Dimensions',{'K' Inf});
    nccreate(obsfile,'type','Dimensions',{'K' Inf},'Datatype','int32');
    nccreate(obsfile,'lon','Dimensions',{'K' Inf});
    nccreate(obsfile,'lat','Dimensions',{'K' Inf});
    nccreate(obsfile,'z','Dimensions',{'K' Inf});
    nccreate(obsfile,'time','Dimensions',{'K' Inf});
    nccreate(obsfile,'dir','Dimensions',{'K',Inf});
    Ktotal=0;

    
    %U,V
    if false
        tLim=tOut(it)+[0,tWin];
        tUV=tOut(it)+[.5:1:tWin];
        obslist_hf(grdfile,fullfile(obsDir,'hfr','kosro/'),obsfile,tLim,tUV,mask);
    end
    if true
        tLim=tOut(it)+[0,tWin]; %time limits
        hfrfile{1}=fullfile(obsDir,'hfr','kosro','RDLi_CBL1.mat');
        hfrfile{2}=fullfile(obsDir,'hfr','kosro','RDLi_LOO1.mat');
        hfrfile{3}=fullfile(obsDir,'hfr','kosro','RDLi_MAN1.mat');
        hfrfile{4}=fullfile(obsDir,'hfr','kosro','RDLi_PSG1.mat');
        hfrfile{5}=fullfile(obsDir,'hfr','kosro','RDLi_WIN1.mat');
        hfrfile{6}=fullfile(obsDir,'hfr','kosro','RDLi_WLD2.mat');
        hfrfile{7}=fullfile(obsDir,'hfr','kosro','RDLi_WSH1.mat');
        hfrfile{8}=fullfile(obsDir,'hfr','kosro','RDLi_YHL1.mat');
        hfrfile{9}=fullfile(obsDir,'hfr','kosro','RDLi_YHS2.mat');
        obs1=obslist_hf_radialsuper(grdfile,hfrfile,tideUV,obsfile,tLim,mask);
    end

    
    %SST
    if true
        tLim=tOut(it)+[0,tWin];
        goesfile={'/home/server/student/homes/ipasmans/Obs/sst/avhrr2011/metopa_l3u_2011.mat',...
            '/home/server/student/homes/ipasmans/Obs/sst/avhrr2011/noaa19_l3u_2011.mat'};
        obs1=obslist_l3u(goesfile,grdfile,obsfile,tLim,...
            mask);
    end
    if false
        tLim=tOut(it)+[0,tWin];
        goesDir='/home/server/student/homes/ipasmans/Obs/sst/goes2011/';
        goesfile={fullfile(goesDir,'goes_l2pf_201104.mat'),...
            fullfile(goesDir,'goes_l2pf_201105.mat'),...
            fullfile(goesDir,'goes_l2pf_201106.mat'),...
            fullfile(goesDir,'goes_l2pf_201107.mat'),...
            fullfile(goesDir,'goes_l2pf_201108.mat'),...
            fullfile(goesDir,'goes_l2pf_201109.mat'),...
            };
        obs1=obslist_goesl2p(goesfile,grdfile,obsfile,tLim,...
            mask);
    end
    if false
         tLim=tOut(it)+[0,tWin]; 
         tGoes=tOut(it)+[0:72]/24; 
         goesDir='/home/server/student/homes/ipasmans/Obs/sst/goes2011/';
         goesfile={fullfile(goesDir,'goesf_20110401_20110601.nc'),...
            fullfile(goesDir,'goesf_20110601_20110801.nc'),...
            fullfile(goesDir,'goesf_20110801_20111001.nc')};
        obslist_goes(goesfile,grdfile,obsfile,tLim,tGoes,mask); 
    end
    if false
        tLim=tOut(it)+[0,tWin];
        %tGoes=tOut(it)+[2/24:4/24:tWin-2/24];
        tGoes=tOut(it)+[0:1/24:tWin];
        goesDir='/home/server/student/homes/ipasmans/Obs/sst/goes2011/';
        goesfile={fullfile(goesDir,'goesf_20110401_20110601.nc'),...
            fullfile(goesDir,'goesf_20110601_20110801.nc'),...
            fullfile(goesDir,'goesf_20110801_20111001.nc')};
        obslist_goes_meanBlock(goesfile,grdfile,obsfile,tLim,60e3,...
                    tGoes,mask);
    end
    if false
        tLim=tOut(it)+[0,tWin];
        goesDir='/home/server/student/homes/ipasmans/Obs/sst/goes2011/';
        goesfile={fullfile(goesDir,'goesl2pf_20110701.mat'),...
            fullfile(goesDir,'goesl2pf_20110801.mat')};
        obslist_goesl2p_meanBlock(goesfile,grdfile,obsfile,tLim,60e3,...
                    mask);
    end
    if false
        tLim=tOut(it)+[0,tWin];
        modisDir='/home/server/student/homes/ipasmans/Obs/sst/modis2011/';
        modisfile={fullfile(modisDir,'modis_20110401_20110601.mat'),...
            fullfile(modisDir,'modis_20110601_20110801.mat'),...
            fullfile(modisDir,'modis_20110801_20111001.mat')};
        obslist_modis(modisfile,grdfile,obsfile,tLim,mask)
    end
    
    %SSH
    if true
        tLim=tOut(it)+[0,tWin];
        obs1=obslist_ssh(grdfile,'/home/server/student/homes/ipasmans/Obs/ssh/',tideZ,obsfile,tLim,mask);
    end
 
    %glider (don't forget to set temperature/salinity in obslist_glider)
    if true
        tLim=tOut(it)+[0,tWin];
        gliderDataDir=fullfile(obsDir,'glider','osu');
        obslist_glider(grdfile,gliderDataDir,obsfile,tLim,mask);
    end

    %check nan
    obs.val=ncread(obsfile,'obs');
    obs.lon=ncread(obsfile,'lon');
    obs.lat=ncread(obsfile,'lat');
    obs.type=ncread(obsfile,'type');
    obs.t=ncread(obsfile,'time');
    obs.sig=ncread(obsfile,'sig_d');
    obs.z=ncread(obsfile,'z'); 
    obs.dir=ncread(obsfile,'dir'); 
    
    %times to remove
    tRemove=[...
    2011 7 21 6 0 0 2011 7 21 13 0 0;...
    2011 7 24 16 0 0 2011 7 24 16 0 0; ...
    2011 7 26 14 0 0 2011 7 26 14 0 0;...
    2011 7 26 16 0 0 2011 7 26 16 0 0;...
    2011 7 27 3 0 0 2011 7 27 3 0 0;...
    2011 7 29 13 0 0 2011 7 29 13 0 0;...
    2011 7 30 4 0 0 2011 7 30 4 0 0;...
    2011 7 30 16 0 0 2011 7 30 16 0 0;...
    2011 7 30 19 0 0 2011 7 30 19 0 0;...
    2011 8 3 3 0 0 2011 8 3 18 0 0;...
    2011 8 5 3 0 0 2011 8 5 17 0 0;...
    2011 8 6 3 0 0 2011 8 6 15 0 0;...
    2011 8 7 14 0 0 2011 8 7 14 0 0;...
    2011 8 8 14 0 0 2011 8 8 18 0 0;...
    2011 8 9 14 0 0 2011 8 9 14 0 0;...
    2011 8 10 6 0 0 2011 8 10 6 0 0;...
    ];  
    tRemove=[...
        2011 08 03 10 0 0 2011 08 03 10 0 0;...
        2011 08 05 15 0 0 2011 08 05 17 0 0;...
        2011 08 06 03 0 0 2011 08 06 10 0 0]; 
        
        
    for itR=[] %1:size(tRemove,1)
        in=obs.type==4 & ...
            obs.t/3600/24+tOut(it)>=datenum(tRemove(itR,1:6)) & ...
            obs.t/3600/24+tOut(it)<=datenum(tRemove(itR,7:end)); 
        obs.lon=obs.lon(~in); obs.lat=obs.lat(~in); obs.z=obs.z(~in); 
        obs.t=obs.t(~in); obs.type=obs.type(~in); obs.sig=obs.sig(~in); obs.val=obs.val(~in); 
    end
       

    if false
        %demean
        for i1=unique(obs.type)
            inan=obs.type==i1;
            obs.val(inan)=obs.val(inan)-mean(obs.val(inan));
        end
    end
    
    
    %check output
    if any(isnan(obs.val))
        error(sprintf('NaN in val %s',obsfile)); 
    end
    if any(isnan(obs.sig))
        error(sprintf('NaN in sig %s',obsfile)); 
    end
    if any(obs.lon>0)
        error(sprintf('Wrong longitude in %s',obsfile)); 
    end
    
    %Write to output
    write_obslist(obsfile,obs); 
   
    display(sprintf('Writing %s DONE',obsfile)); 
   
end


